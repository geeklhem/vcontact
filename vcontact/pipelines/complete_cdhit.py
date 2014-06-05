"""
STEPS :
0. pblast : Run a pblast all prots against all prots (with refesq as defined in options.refseq.faa (saved in pblast/)
1. mcl : Run mcl to build the protein clusters. (saved in protein_clusters/)
2. network : Load mcl results into python, build the hypergeometric contig network (saved in networks/) 
3. affiliate : Run mcl to build contigs clusters. Load results into python 
"""

import logging
import os
import re 
import argparse
import cPickle as pickle
import subprocess
import vcontact.options as options
import vcontact.genomes as genomes
import vcontact.protein_clusters as protein_clusters
import vcontact.genome_clusters as genome_clusters
import vcontact.pc_matrix as pc_matrix


step_list = {"pblast":0,
             "mcl":1,
             "network":2,
             "affiliate":3,
             "nostop":100}


# Argparse config 
parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter, )
parser.add_argument('step',
                    help='Starting step of the pipeline',
                    choices=step_list.keys(),
                    nargs="?",
                    default=[k for k, v in sorted(step_list.items(), key=lambda x:x[1])][0])
parser.add_argument('-i', '--input',
                    nargs="*",
                    help='protein sequence files', metavar='FAA')
parser.add_argument('-f', '--force-overwrite',
                    help='Overwrite existing files',
                    action="store_true")
parser.add_argument('--tara',
                    help='Load tara contigs and metadata',
                    action="store_true")
parser.add_argument('--refcont',
                    help='Load refseq simulated contigs metadata',
                    action="store_true")
parser.add_argument('-v', '--verbose',
                    help='Verbosity level : -v warning, -vv info, -vvv debug, (default debug)',
                    action="count",
                    default=-1)
parser.add_argument('-c', '--identity',
                    help='CDHIT Identity threshold',
                    default=0.9,
                    type=float)
parser.add_argument('-i2','--inflation2',
                    help='MCL Inflation',
                    default=2)
parser.add_argument('--stop',
                    help='Stopping step',
                    default="nostop")


args = parser.parse_args()


# Logging config
log_levels = [logging.WARNING, logging.INFO, logging.DEBUG]

# create logger with 'spam_application'
logger = logging.getLogger('vcontact')
logger.setLevel(logging.DEBUG)
# create file handler which logs even debug messages
fh = logging.FileHandler(options.folders["log"] + "complete_pipeline.log")
fh.setLevel(logging.INFO)
# create console handler with a higher log level
ch = logging.StreamHandler()
ch.setLevel(log_levels[args.verbose])
# create formatter and add it to the handlers
fh.setFormatter(logging.Formatter('%(asctime)s - %(process)d - %(levelname)s - %(message)s'))
ch.setFormatter(logging.Formatter('%(levelname)s - %(message)s'))
# add the handlers to the logger
logger.addHandler(fh)
logger.addHandler(ch)

logger.info("{:#^80}".format(""))
logger.info("{:.^80}".format(" AUTOMATIC TAXONOMY PIPELINE  "))
logger.info("Started at step {} @ {}".format(args.step, " ".join(os.uname())))


g = genomes.Genomes(refseq=True, tara=args.tara, refseq_contigs=args.refcont)


def call(cmd, fi, overwrite=None):
    """
    Call the command 'cmd' (using subprocess) if the file 'fi' does not exist.
    The file 'fi' is removed if the command 'cmd' fails or is interrupted by the user.

    - 'cmd' (str) : shell cmd. 
    - 'fi' (str) : file (should be the output file of 'cmd').
    - 'overwrite' (bool) : If True, 'cmd' is executed even if 'fi' exist.

    Security concern : This method is NOT bash-injection safe. 
    """
    overwrite = options.overwrite if overwrite == None else overwrite
    if overwrite or not os.path.isfile(fi):
        logger.info("Executing {} to create file {}.".format(cmd.split(" ")[0], os.path.basename(fi)))
        try:
            subprocess.check_call(cmd, shell=True) #Security issue : vulnerable to shell injections...
        except subprocess.CalledProcessError as e:
            logger.error("Command {} failed ({}). Removing {}.".format(cmd.split(" ")[0], e, os.path.basename(fi)))
            os.remove(fi)
            raise
        except KeyboardInterrupt as e:
            logger.error("Command {} interrupted. Removing {}.".format(cmd.split(" ")[0], os.path.basename(fi)))
            os.remove(fi)
            raise 
            
    else:
        logger.info("File {} exist, command {} aborted.".format(fi, cmd.split(" ")[0]))

    if not os.path.isfile(fi):
        logger.info("{1} failed to create file {}.".format(fi, cmd.split(" ")[0]))
        raise IOError


def wordsize(threshold):
    """Given the identity threshold, return the n parameter value for cdhit. """
    if 0.7 <= threshold <= 1:
        return 5
    elif 0.6 <= threshold < 0.7:
        return 4
    elif 0.5 <=threshold < 0.6:
        return 3
    elif 0.4 <= threshold < 0.5:
        return 2
    else:
        raise ValueError("Invalid threshold".format(threshold))


if __name__ == "__main__":
    
    os.chdir(options.data_folder)
    
    for fi in args.input:
        name = "_".join(os.path.basename(os.path.expanduser(fi)).split(".")[:-1]) 
        path = options.folders["results"]+name+"/"
        if not os.path.exists(path):
            os.mkdir(path)
            logging.info("Created {}/ directory into {}".format(name, options.folders["results"]))

        logger.info("{:.^80}".format("file : {} ".format(fi)))

        ### Files ###
        # - fi, refseq_faa
        # - complete_faa  : concatenation of fi and refseq_faa
        # - output_cdhit  : files outputed by cdhit (two of them .faa & .faa.clstr)
        # - protein_clusters : a line by protein cluster

        cdhit_param = "cdhit{}".format("".join(str(args.identity).split(".")))
        param = "i{}_cdhit{}".format("".join(str(args.inflation2).split(".")),
                                     "".join(str(args.identity).split(".")))
        files = {
            'input': os.path.expanduser(fi),
            'refseq_faa': options.refseq_faa,
            'proteins_faa' : path+"proteins.faa",
            'output_cdhit' : path+"proteins_{}".format(cdhit_param),
            'protein_clusters' : path+"proteins_{}.clusters".format(cdhit_param), 
            "pc_matrix" : path+"pcm_{}.pkle".format(cdhit_param),
            "ntw": path+"contigs_{}.ntw".format(cdhit_param),
            "gc": path+"gc_{}.pkle".format(param),
            "contig_clusters":path+"contigs_{}.clusters".format(param),
            "contigs":path+"contigs_{}.pandas".format(param),
            "clusters":path+"clusters_{}.pandas".format(param),
            "summary":path+"summary_{}.pandas".format(param),
        }


        ### STEP 1: CLUSTERING 
        if step_list[args.step]<2<=step_list[args.stop]+1:
            print("{0:-^80}".format(" CDHIT "))

            ## Concatenate 
            logger.info("Concatenation proteins fasta files...")
            call("awk {{print}} {refseq} {i} > {o}".format(o=files["proteins_faa"],
                                                           i=files["input"],
                                                           refseq=files["refseq_faa"]),
                 files["proteins_faa"])

            ## CD-hit
            logger.info("Running CDhit...")
            call("cd-hit -i {i} -o {o}.faa -c {c} -n {n} -d 0".format(i=files["proteins_faa"],
                                                                      o=files["output_cdhit"],
                                                                      c=args.identity,
                                                                      n=wordsize(args.identity)),
                 files["output_cdhit"]+".faa")


            ## Parsing CD-hit results
            logger.info("Parsing CD-hit results...")
            l = 0
            with open(files["output_cdhit"]+".faa.clstr", 'r') as fi_clstr:
                with open(files["protein_clusters"], "w") as fi_clusters:
                    for line in fi_clstr:
                        if re.match("^>Cluster", line):
                            fi_clusters.write("\n")
                            l = l+1
                        else:
                            fi_clusters.write(line.split(">")[1].split("...")[0])
            logger.info("CDHIT completed : {} proteins clusters created.".format(l))


        ### STEP 2: NETWORK
        if step_list[args.step]<3<=step_list[args.stop]+1:
            print("{0:-^80}".format(" NETWORK "))

            ## Load the pcs
            logger.info("Loading the protein clusters...")
            p = protein_clusters.ProteinClusters(files["protein_clusters"])



            ## Building the network 
            if options.overwrite or not os.path.isfile(files["pc_matrix"]):
                logger.info("Building the network")
                pcm = pc_matrix.PCMatrix(p.data.proteins, p.data.clusters,
                                         g.data.proteins, g.data.contigs)

                with open(files["pc_matrix"], 'wb') as f:
                    pickle.dump(pcm, f)
            else:
                logger.info("Loading the network")
                with open(files["pc_matrix"], 'rb') as f:
                    pcm = pickle.load(f)

            ## Exporting the network 
            logger.info("Exporting the network")
            if options.overwrite or not os.path.isfile(files["ntw"]):
                pcm.to_mcl(pcm.ntw, files["ntw"])

        ### STEP 3: CONTIGS CLUSTERS
        if step_list[args.step]<4<=step_list[args.stop]+1:
            print("{0:-^80}".format(" CONTIGS CLUSTERS "))


            ## Clustering the contigs
            logger.info("Clustering the contigs")
            call("mcl {0} --abc -I {1} -o {2} ".format(files["ntw"], args.inflation2, files["contig_clusters"]),
                 files["contig_clusters"])

            ## Analysis of the network
            if options.overwrite or not os.path.isfile(files["gc"]):
                logger.info("Computing the contigs clusters and affiliation")
                gc = genome_clusters.GenomeCluster(pcm,
                                                   mcl_file=files["contigs_clusters"],
                                                   name="".join(os.path.basename(files["contigs_clusters"]).split(".")[:-1]))
                gc.routine()
                gc.nodes_properties()
                gc.nodes_size()

                # Exports
                logger.info("Exporting")
                gc.contigs.to_pickle(files["contigs"])
                gc.clusters.to_pickle(files["clusters"])
                gc.summary.to_pickle(files["summary"])
                with open(files["gc"], "w") as f:
                    pickle.dump(files["gc"], f)

            else:
                logger.info("Loading the contigs clusters properties")
                with open(files["gc"], 'rb') as f:
                    gc = pickle.load(f)



