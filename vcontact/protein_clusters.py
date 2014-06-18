"""Protein_clusters.py"""

import pandas
import logging
import os
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
        print [p for p in prots if p not in proteins.index] 
        prots = [p for p in prots if p in proteins.index] 
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

def load_refseq(self,fi=None,):
        """ Load the data from the genbank file to creat the protein.csv, contigs.csv files.

        Args:
            fi (str): The genebank file.
            fi_taxt (str): a pickled pandas dataframe with the taxonomy information.
            h5_name (str): path to the cache file to save. 
        
        Returns:
           (str) proteins.csv, contigs.csv 
        """

        # Loading taxonomic classes 
        fi = os.path.expanduser(options.data_folder+"refseq/phage.genomic.gbff") if not fi else fi 
        fi_taxa = os.path.expanduser(options.data_folder+"ncbi_taxdump/refseqphage_taxa.pkl") if not fi_taxa else fi_taxa  


        store =  pandas.HDFStore(options.cache_folder+h5_name)

        if "contigs" not in store or "proteins" not in store:
            taxa = pandas.read_pickle(fi_taxa)
            families = frozenset(taxa[taxa["Rank"]=="family"]["Name"])
            genera = frozenset(taxa[taxa["Rank"]=="genus"]["Name"])


            genome = {"name":[],
                      "family":[],
                      "genus":[],
                      "size":[],}

            proteins = {"protein_id":[],
                        "contig":[],
                        "function":[],}

            # Reading the gbff file thanks to BioPython
            rec = SeqIO.parse(fi, "gb")
            for r in rec:
                name = r.id.split(".")[0]
                genome["name"].append(name)

                # Taxonomy
                genome["size"].append(len(r.seq))
                genome["family"].append(None)
                genome["genus"].append(None)
                for t in r.annotations["taxonomy"]:
                    if t in families:
                        genome["family"][-1] = t
                    elif t in genera:
                        genome["genus"][-1]  = t
                #Proteins 
                for feature in r.features:
                    if "protein_id" in feature.qualifiers:
                        proteins["protein_id"].append(feature.qualifiers["protein_id"][0].split(".")[0])
                        proteins["contig"].append(name)
                        proteins["function"].append(feature.qualifiers["product"][0])

                 

            genome = pandas.DataFrame(genome).set_index("name")
            genome["origin"] = "refseq_jan14"
            proteins = pandas.DataFrame(proteins).set_index("protein_id")

            pbc = pandas.DataFrame(proteins.groupby("contig").contig.count(), columns=["proteins"])
            genome = pandas.merge(genome,pbc,how="left",left_index=True,right_index=True)


            store.append('contigs',genome,format='table',data_columns=True)
            store.append('proteins',proteins,format='table',data_columns=True)

            nf = len(genome) - genome.family.count() 
            ng = len(genome) - genome.genus.count()
            if nf!=0 or ng!=0 :
                logger.warning("{0} genomes without family and {1} genomes without genus.".format(nf,ng))

        return store 

