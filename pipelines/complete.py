"""
STEPS :
0. pblast : Run a pblast all prots against all prots (with refesq as defined in options.refseq.faa (saved in pblast/)
1. mcl : Run mcl to build the protein clusters. (saved in protein_clusters/)
2. network : Load mcl results into python, build the hypergeometric contig network (saved in networks/) 
3. affiliate : Run mcl to build contigs clusters. Load results into python 
"""

import logging
import platform, sys, os
import glob, re 
import argparse
import cPickle as pickle
import subprocess
import vcontact.options as options
import vcontact.genomes as genomes
import vcontact.protein_clusters as protein_clusters
import vcontact.genome_clusters as genome_clusters
import vcontact.pc_matrix as pc_matrix
import time

step_list = {"pblast":0,
             "mcl":1,
             "network":2,
             "affiliate":3}


# Argparse config 
parser = argparse.ArgumentParser(description=__doc__,formatter_class=argparse.RawDescriptionHelpFormatter,)
parser.add_argument('step',
                    help='Starting step of the pipeline',
                    choices=step_list.keys(),
                    nargs="?",
                    default=[k for k,v in sorted(step_list.items(),key=lambda x:x[1])][0])
parser.add_argument('-i','--input',
                    nargs="*",
                    help='protein sequence files',metavar='FAA')
parser.add_argument('-f','--force-overwrite',
                    help='Overwrite existing files',
                    action="store_true")
parser.add_argument('--tara',
                    help='Load tara contigs and metadata',
                    action="store_true")
parser.add_argument('-v','--verbose',
                    help='Verbosity level : -v warning, -vv info, -vvv debug, (default debug)',
                    action="count",
                    default=-1)
parser.add_argument('-i2','--inflation2',
                    help='Inflation of the second MCL',
                    default=2)


args = parser.parse_args()

# Logging config
log_levels = [logging.WARNING,logging.INFO,logging.DEBUG]

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
logger.info("Started at step {} @ {}".format(args.step," ".join(os.uname())))


g = genomes.Genomes(refseq=True,tara=args.tara)

def call(cmd,fi,overwrite=None):
    overwrite = options.overwrite if overwrite == None else overwrite
    if overwrite or not os.path.isfile(fi):
        logger.info("Executing {} to create file {}.".format(cmd.split(" ")[0],fi))
        subprocess.check_call(cmd, shell=True) #Security issue : vulnerable to shell injections...
    else:
        logger.info("File {} exist, command {} aborted.".format(fi,cmd.split(" ")[0]))

if __name__ == "__main__":
    
    os.chdir(options.data_folder)
    for fi in args.input:
        try:
            name = "_".join(os.path.basename(fi).split(".")[:-1])+"_and_refseq"
            logger.info("{:.^80}".format("file : {} ".format(fi)))

            ### STEP 0 : PBLAST 
            if step_list[args.step]<1:
                print("{0:-^80}".format(" AvA PBLAST "))

                logger.info("  Concatenate with refseq")
                faa = options.folders["proteins"]+name+".faa"
                call("cat {i} {r} > {o}".format(i=fi,r=options.refseq_faa,o=faa),faa)

                logger.info("  Creating a database")
                db = options.folders["proteins"]+name+".db"
                call("formatdb -i {i} -n {o}".format(i=faa,o=db),db+".phr")

                logger.info("  All versus all blast")
                blasted = options.folders["proteins"]+name+".tab"
                call("blastall -p blastp -i {i} -d {db} -o {o} -m 8 -a 24 -e 0.00001".format(i=faa,
                                                                                             db=db,
                                                                                             o=blasted),
                     blasted)

            ### STEP 1: MCL
            if step_list[args.step]<2:
                logger.info("  MCL")
                path = options.folders["proteins"]+name
                call("awk '$1!=$2 {{print $1,$2,$11}}' {0}.tab > {0}.abc".format(path),path+".abc")
                call("mcxload -abc  {0}.abc --stream-mirror --stream-neg-log10 -stream-tf 'ceil(200)' \
                -o {0}.mci -write-tab {0}_mcxload.tab".format(path),path+".mci")
                call("mcl {0}.mci -I 2 -use-tab {0}_mcxload.tab -o {0}.clusters".format(path),path+".clusters")

            ### STEP 2: NETWORK
            if step_list[args.step]<3:
                logger.info("  Network")
                p = protein_clusters.ProteinClusters(options.folders["proteins"]+name+".clusters")
                path = options.folders["contigs"]+name
                pcm_pklefile = path+"_pc_matrix.pkle"
                if options.overwrite or not os.path.isfile(pcm_pklefile):
                    pcm = pc_matrix.PCMatrix(p.data.proteins,p.data.clusters,
                                             g.data.proteins, g.data.contigs)

                    with open(pcm_pklefile, 'wb') as f:
                        pickle.dump(pcm,f)
                else:
                    with open(pcm_pklefile, 'rb') as f:
                        pcm = pickle.load(f)

                if options.overwrite or not os.path.isfile(path+".ntwk"):
                    pcm.to_mcl(pcm.ntw,"{0}.ntwk".format(path))

            ### STEP 3: CONTIGS CLUSTERS
            if step_list[args.step]<4:
                path = options.folders["contigs"]+name
                clustersfile = {0}_i_{0}_contigs.clustersformat(path,"".join(str(args.inflation2).split(".")))
                call("mcl {0}.ntwk --abc -I {1} -o {2} ".format(path,,args.inflation2,clustersfile),clustersfile)
                gc_pklefile = path+"_contig_clusters.pkle"
                if options.overwrite or not os.path.isfile(gc_pklefile):
                    gc = genome_clusters.GenomeCluster(pcm,mcl_file=path+"_contigs.clusters")
                    gc.routine()
                    with open(gc_pklefile, 'wb') as f:
                        pickle.dump(gc,f)
                else:
                    with open(gc_pklefile, 'rb') as f:
                        gc = pickle.load(f)
        except Exception:
            logging.error("Error in {}".format(fi))



