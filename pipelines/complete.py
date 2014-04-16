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
             "affiliate":3,
             "nostop":100}


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
parser.add_argument('--refcont',
                    help='Load refseq simulated contigs metadata',
                    action="store_true")
parser.add_argument('-v','--verbose',
                    help='Verbosity level : -v warning, -vv info, -vvv debug, (default debug)',
                    action="count",
                    default=-1)
parser.add_argument('-i2','--inflation2',
                    help='Inflation of the second MCL',
                    default=2)
parser.add_argument('--stop',
                    help='Stopping step',
                    default="nostop")


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


g = genomes.Genomes(refseq=True,tara=args.tara,refseq_contigs=args.refcont)

def call(cmd,fi,overwrite=None):
    overwrite = options.overwrite if overwrite == None else overwrite
    if overwrite or not os.path.isfile(fi):
        logger.info("Executing {} to create file {}.".format(cmd.split(" ")[0],os.path.basename(fi)))
        try:
            subprocess.check_call(cmd, shell=True) #Security issue : vulnerable to shell injections...
        except subprocess.CalledProcessError as e:
            logger.error("Command {} failed ({}). Removing {}.".format(cmd.split(" ")[0],e,os.path.basename(fi)))
            os.remove(fi)
            raise
        except KeyboardInterrupt as e:
            logger.error("Command {} interrupted. Removing {}.".format(cmd.split(" ")[0],os.path.basename(fi)))
            os.remove(fi)
            raise 
            
    else:
        logger.info("File {} exist, command {} aborted.".format(fi,cmd.split(" ")[0]))

    if not os.path.isfile(fi):
        logger.info("{1} failed to create file {}.".format(fi,cmd.split(" ")[0]))
        raise IOError

if __name__ == "__main__":
    
    os.chdir(options.data_folder)
    refseq_faa = options.refseq_faa
    print("{0:-^80}".format(" Refseq  "))
    
    logger.info("  Creating a database")
    refseq_db = options.folders["proteins"]+"refseq.db"
    call("formatdb -i {i} -n {o}".format(i=refseq_faa,o=refseq_db),refseq_db+".phr")
    
    logger.info("  All versus all blast")
    refseq_blasted = options.folders["proteins"]+"refseq.tab"
    call("blastall -p blastp -i {i} -d {db} -o {o} -m 8 -a 24 -e 0.00001".format(i=refseq_faa,
                                                                                 db=refseq_db,
                                                                                 o=refseq_blasted),
         refseq_blasted)


    
    for fi in args.input:
        try:
            fi = os.path.expanduser(fi)


            faa_basename = "_".join(os.path.basename(fi).split(".")[:-1])
            name = "_".join(os.path.basename(fi).split(".")[:-1])+"_and_refseq"
            
            logger.info("{:.^80}".format("file : {} ".format(fi)))

            ### STEP 0 : PBLAST 
            if step_list[args.step]<1<=step_list[args.stop]+1:
                print("{0:-^80}".format(" PBLAST "))

                ### STEP 0.1 : CONTIGS VS REF SEQ
                logger.info("  Contigs VS Refseq")
                contigs_v_refseq = options.folders["proteins"]+faa_basename+"_v_refseq.tab"
                
                call("{cmd} -i {i} -d {db} -o {o}".format(cmd=options.blast_cmd,
                                                          i=fi,
                                                          db=refseq_db,
                                                          o=contigs_v_refseq),
                     contigs_v_refseq)

                ### STEP 0.2 : CONTIGS + REF SEQ VS CONTIGS
                logger.info("  Contigs+Refseq VS Contigs")
                c_and_r_v_c = options.folders["proteins"]+faa_basename+"_and_refseq_v_contg.tab"
                db = options.folders["proteins"]+faa_basename+".db"
                merged_fi = options.folders["proteins"]+faa_basename+"_and_refseq.faa"
                
                call("formatdb -i {i} -n {o}".format(i=fi,o=db),db+".phr")
                call("cat {i} {r} > {o}".format(i=fi,r=options.refseq_faa,o=merged_fi),merged_fi)
                call("{cmd} -i {i} -d {db} -o {o}".format(cmd=options.blast_cmd,
                                                          i=merged_fi,
                                                          db=db,
                                                          o=c_and_r_v_c),
                     c_and_r_v_c)
                
                ### STEP 0.3 : MERGING THE RESULTS
                logger.info("  Merging results RvR, CvR and (C+R)vC")
                blasted = options.folders["proteins"]+name+".tab"
                blast_files = [refseq_blasted, contigs_v_refseq, c_and_r_v_c]
                call("cat {i} > {o}".format(i=" ".join(blast_files),o=blasted), blasted)

            ### STEP 1: MCL
            if step_list[args.step]<2<=step_list[args.stop]+1:
                print("{0:-^80}".format(" MCL "))
                path = options.folders["proteins"]+name
                call("awk '$1!=$2 {{print $1,$2,$11}}' {0}.tab > {0}.abc".format(path),path+".abc")
                call("mcxload -abc  {0}.abc --stream-mirror --stream-neg-log10 -stream-tf 'ceil(200)' \
                -o {0}.mci -write-tab {0}_mcxload.tab".format(path),path+".mci")
                call("mcl {0}.mci -I 2 -use-tab {0}_mcxload.tab -o {0}.clusters".format(path),path+".clusters")

            ### STEP 2: NETWORK
            if step_list[args.step]<3<=step_list[args.stop]+1:
                print("{0:-^80}".format(" NETWORK "))
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
            if step_list[args.step]<4<=step_list[args.stop]+1:
                print("{0:-^80}".format(" CONTIGS CLUSTERS "))
                path = options.folders["contigs"]+name
                clustersfile = "{0}_i_{1}_contigs.clusters".format(path,"".join(str(args.inflation2).split(".")))
                call("mcl {0}.ntwk --abc -I {1} -o {2} ".format(path,args.inflation2,clustersfile),clustersfile)
                gc_pklefile = path+"_contig_clusters.pkle"
                if options.overwrite or not os.path.isfile(gc_pklefile):
                    gc = genome_clusters.GenomeCluster(pcm,
                                                       mcl_file=clustersfile,
                                                       name="".join(os.path.basename(clustersfile).split(".")[:-1]))
                    gc.routine()
                    with open(gc_pklefile, 'wb') as f:
                        pickle.dump(gc,f)
                else:
                    with open(gc_pklefile, 'rb') as f:
                        gc = pickle.load(f)
        except Exception as e:
            logger.error("Error {} in {}".format(e,os.path.basename(fi)))



