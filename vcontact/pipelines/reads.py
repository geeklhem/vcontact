# coding: utf-8
import logging 
import cPickle as pickle
import glob
from Bio import SeqIO
import numpy as np
import os 

import options
import genome_clusters 
import pc_matrix 
import pandas
import subprocess
import exports

logging.basicConfig(format='%(levelname)s: %(message)s', level=logging.DEBUG)

folder = options.data_folder + "reads/"
fi_size = options.data_folder + "cache/tara_contigs_size.pandas"
fi_size_csv = options.data_folder + "cache/tara_contigs_size.csv"

if not os.path.exists(fi_size):
   
    if not os.path.exists(fi_size_csv):
        logging.info("Computing the contig sizes...")
        rec = SeqIO.parse(options.data_folder+"tara/ALL_CONTIGS_NO_QC.FNA","fasta")
        with open(fi_size_csv,"w") as f:
            f.write("name\tsize\n")
            for i,r in enumerate(rec):
                if i%100000==0:
                    logging.debug("{:0.2%}".format(i/3821756.0))
                f.write("\t".join([str(x) for x in [r.id,len(r.seq)]]))
                f.write("\n")

    logging.info("Origin extraction...")
    df_size = pandas.read_csv(fi_size_csv,sep="\t")
    df_size["origin"] = [x.split("_")[0] for x in df_size.name]
    df_size["station"] = [x[:-3] for x in df_size.origin]
    df_size["zone"] = [x[-3:] for x in df_size.origin]
    logging.info("Saving...")
    df_size.to_pickle(fi_size)

else:
    logging.info("Loading the contig sizes")
    df_size = pandas.read_pickle(fi_size)

for i,fi in enumerate(glob.glob(folder+"*.TAB")):
    logging.info("{:0.2%} - {}".format(i/43.0,fi))
    fi_out = "".join(fi.split(".")[:-1])

    a = pandas.read_csv(fi, sep="\t", header=None, names=["name","reads"])
    a = pandas.merge(a, df_size, how="left")

    a["abundance"] = a.reads/a.reads.sum()
    a.abundance = a.reads/a.loc[:,"size"]
    a["proportions"] = a.abundance/a.abundance.sum()

    a.to_csv(fi_out+".csv")
    a.to_pickle(fi_out+".pandas")
