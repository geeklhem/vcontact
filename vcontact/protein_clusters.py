"""Protein_clusters.py"""

import pandas
import logging
import os
import options
import re 
import subprocess
import scipy.sparse as sparse 
import numpy as np
import sys
import cPickle as pickle

logger = logging.getLogger(__name__)


def make_protein_clusters(blast_fi,path,inflation=2):
    """
    Args: 
        blast_fi (str): Blast results file
        inflation (float): MCL's inflation
        path (str): file basename.
    Returns:
        str: MCL clustering file.  
    """

    logger.debug("Generating abc file...")
    subprocess.check_call(("awk '$1!=$2 {{print $1,$2,$11}}' {0} "
                           "> {1}.abc").format(blast_fi,path),shell=True)
    logger.debug("Running MCL...")
    subprocess.check_call(("mcxload -abc  {0}.abc --stream-mirror "
                           "--stream-neg-log10 -stream-tf 'ceil(200)' " 
                           "-o {0}.mci -write-tab {0}_mcxload.tab "
                           "").format(path),shell=True)
    subprocess.check_call(("mcl {0}.mci -I {1} -use-tab "
                           "{0}_mcxload.tab -o {0}.clusters").format(path,inflation),shell=True)

    return path+".clusters"


def load_clusters(fi,proteins):
    """
    Load given clusters file
    
    Args:
        fi (str): path to clusters file
        proteins (dataframe): A dataframe giving the protein and its contig.
    Returns: 
        tuple: dataframe proteins and dataframe clusters
    """
    
    # Read MCL
    with open(fi) as f:
            c = [ line.rstrip("\n").split("\t") for line in f ]
    c = [x for x in c if len(c)>1]
    nb_clusters = len(c)
    formater = "PC_{{:>0{}}}".format(int(round(np.log10(nb_clusters))+1))
    name = [formater.format(str(i)) for i in range(nb_clusters)]
    size = [len(i) for i in c]
    clusters = pandas.DataFrame({"size":size,"id":name}).set_index("id")
    
    # Assign each prot to its cluster 
    proteins = proteins.set_index("id")
    for prots,clust in zip(c,name):
        proteins.loc[prots,"cluster"] = clust 

    # Keys
    for clust,prots in proteins.groupby("cluster"):
        clusters.loc[clust,"annotated"] = prots.keywords.count()
        if prots.keywords.count():
            keys = ";".join(prots.keywords.dropna().values).split(";")
            key_count = {}
            for k in keys:
                k = k.strip()
                try:
                    key_count[k] += 1
                except KeyError:
                    key_count[k] = 1
            clusters.loc[clust,"keys"] = "; ".join(["{} ({})".format(x,y) for x,y in key_count.items()])

    proteins.reset_index(inplace=True)
    clusters.reset_index(inplace=True)
    profiles = proteins.loc[:,["contig","cluster"]].drop_duplicates()
    profiles.columns = ["contig_id","pc_id"]
    contigs = pandas.DataFrame(proteins.fillna(0).groupby("contig").count().contig)
    contigs.index.name = "id"
    contigs.columns=["proteins"]
    contigs.reset_index(inplace=True)
    
    return proteins,clusters,profiles,contigs

