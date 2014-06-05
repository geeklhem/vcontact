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

# logging.basicConfig(format='%(levelname)s: %(message)s', level=logging.DEBUG)

# print "!"

# logging.info("Load affiliation")
# folder=options.data_folder+"tara/reads/"
# aff = pandas.read_pickle(folder+"affiliation_sig1.pandas")
# folder = options.data_folder + "reads/"

def abundance_by_category(aff,col_category="predicted_family",folder=options.data_folder+"reads/"):
    """
    Get the normalised abundance by family of contig.
    
    Args:
        aff (dataframe): with column name and category. 
        col_category (str): category to sort contigs (in aff)
        folder (str): A folder where there is a pandas file by station. 
            With columns name and abundance

    Returns:
        dataframe: With a line by station an a column by category. 
    """
    aff = aff.loc[:,["name",col_category]]
    s = []
    for i,fi in enumerate(glob.glob(folder+"*.pandas")):
        logging.info("{:0.2%} - {}".format(i/43.0,fi))
        a = pandas.read_pickle(fi)
        b = pandas.merge(a, aff, left_on="name",
                         right_on="name", how="left")
        line = b.groupby(col_category).sum().abundance
        line.name = os.path.basename(fi).split(".")[0]
        s.append(line)
    ss = pandas.DataFrame(s)
    return ss

#+ TODO : REFACTOR AFTER THIS LINE
def ab_by_name(dict_names,aff):
    output = dict([(n,[]) for n in dict_names.keys()])
    output["name"] = []
    output["others"] = [] 

    for i,fi in enumerate(glob.glob(folder+"*.pandas")):
        logging.info("{:0.2%} - {}".format(i/43.0,fi))
        output["name"].append(os.path.basename(fi).split(".")[0])
        a = pandas.read_pickle(fi)
        b = pandas.merge(a, aff, left_on="name", right_on="contig", how="left")
        b = b.dropna(axis=0, subset=["family"])
        
        b["c10_proportions"] = b.abundance/b.abundance.sum()

        s= 0
        for k,v in dict_names.items():
            p = np.nan_to_num(float(b.set_index("name").ix[v,"c10_proportions"].sum()))
            output[k].append(p)
            s += p 
        output["others"].append(1-s)
    output = pandas.DataFrame(output)
    return output




