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

folder = options.data_folder + "tara/network/"
contigs = pandas.read_pickle(folder+"contigs_sig1.pandas")
cluster = pandas.read_pickle(folder+"cluster_sig1.pandas")
taxonomy = pandas.read_pickle(folder+"taxonomy_sig1.pandas")

m1 = pandas.merge(contigs.reset_index(),
                  cluster,
                  left_on="pos_cluster",right_on="pos",
                  how="left",
                  suffixes=["__contig","__cluster"])

m2 = pandas.merge(m1,
                 taxonomy.query("level=='family'").loc[:,["pos","name"]],
                  left_on="pos_family", right_on="pos",
                  how="left",
                  suffixes=["","__family"])

aff = m2.loc[:,["name__contig","name"]]
aff.columns = ["contig","family"]
aff.family = aff.family.fillna("Non affiliated")
print aff.groupby("family").count()
aff.to_pickle(folder+"affiliation_sig1.pandas")
