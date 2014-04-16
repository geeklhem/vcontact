"""Protein_clusters.py"""

import pandas
import logging
import os
import options
import re 
from genomes import pid as pid  

class ProteinClusters(object):
    """
    A class to store protein clusters
    """
    
    def __init__(self,clusters_fi):
        """
        """
        self.data = self.load_clusters(clusters_fi)

    def __len__(self):
        return len(self.data.clusters)

    def load_clusters(self,fi,):
        """
        Load given clusters file
        Arguments:
        - `fi`: clusters file
        OUTPUT :
        - (pandas hdf5 store) containing the dataframe proetin and clusters.
        """
        h5_name = ''.join(os.path.basename(fi).split(".")[:-1])+".h5"
        store =  pandas.HDFStore(options.cache_folder+h5_name)

        if "proteins" not in store or "clusters" not in store:  
            proteins = {"protein_id":[],
                        "cluster":[]}
            clusters = {"size":[],
                        "name":[]}
            
            with open(fi) as f:
                for C,l in enumerate(f):
                    l = l.split()
                    clusters["name"].append("pc_{}".format(C))
                    clusters["size"].append(len(l))
                    proteins["cluster"] += ["pc_{}".format(C)] * len(l)
                    try :
                        proteins["protein_id"] += [pid(prot) for prot in l]
                        
                    except AttributeError as e:
                        print prot
                        raise
                                        
            store.append("proteins",pandas.DataFrame(proteins).set_index("protein_id"), format="table")
            store.append("clusters",pandas.DataFrame(clusters).set_index("name"), format="table")

        return store 

    
        
