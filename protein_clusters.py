"""Protein_clusters.py"""

import pandas
import logging
import os
import options
import re 
from genomes import pid as pid  
import genomes
import subprocess
import scipy.sparse as sparse 
import numpy as np
import sys
import cPickle as pickle

logger = logging.getLogger(__name__)

class ProteinClusters(object):
    """
    A class to store protein clusters
    """
    
    def __init__(self,clusters_fi,keys=False):
        """
        Initialise the onject with a cluster file.

        Args:
            cluster_fi (str): path to clusters file
            keys (bool): parse the functions of the proteins in each cluster.
        """
        self.data,self.key_matrix = self.load_clusters(clusters_fi,keys)

    def __len__(self):
        return len(self.data.clusters)

    def load_clusters(self,fi,keys):
        """
        Load given clusters file
        
        Args:
            fi (str): path to clusters file
            keys (bool): parse the functions of the proteins in each cluster
        
        Returns: 
            pandas.HDF5Store: containing the dataframe proetin and clusters.

        Warning:
            The function is painfully slow with keys=True. I'm not sure why and I
            do not have the time to optimize it. It is probably easy to fix though.
        """
        
        g = genomes.Genomes(True,False)
        dataframe_prot = g.data.proteins.copy()
        
        dataframe_prot["anotated"] = [True   if (not len(func) == 0
                                                   and not re.match("[0-9]*ORF[0-9]+",func,flags=re.IGNORECASE) 
                                                   and not re.match("gp[0-9]+",func,flags=re.IGNORECASE) 
                                                   and not re.match("hypothetical",func,flags=re.IGNORECASE))
                                        else False
                                        for func
                                        in dataframe_prot.function]
        dataframe_prot = dataframe_prot.query("anotated==True")
        lenanot = float(len(dataframe_prot))
        
        print "{} annotated proteins".format(lenanot)

        h5_name = ''.join(os.path.basename(fi).split(".")[:-1])+".h5"
        keyfile = ''.join(os.path.basename(fi).split(".")[:-1])+"_keywords.pkle"
        store =  pandas.HDFStore(options.cache_folder+h5_name)
        keywords = pandas.DataFrame({"pos":range(len(options.keywords)),"keyword":options.keywords})
        key_count = np.zeros(len(options.keywords))
        queries = zip(keywords.pos,keywords.keyword)
        
        if "proteins" not in store or "clusters" not in store:  
            proteins = {"protein_id":[],
                        "cluster":[]}
            clusters = {"size":[],
                        "name":[],
                        "annotated":[]}
            
            jj,jjj = 0,0

            with open(fi) as f:
                nb_clusters = len(f.readlines())
                f.seek(0)
                key_matrix = sparse.lil_matrix((len(keywords),nb_clusters))

                
                for C,l in enumerate(f):
                    l = l.split()
                    if len(l) > 1: # drop the singletons
                        clusters["name"].append("pc_{}".format(C))
                        clusters["size"].append(len(l))
                        clusters["annotated"].append(0)
                        proteins["cluster"] += ["pc_{}".format(C)] * len(l)
                        try :
                            proteins["protein_id"] += [pid(prot) for prot in l]
                        except AttributeError as e:
                            print prot
                            raise

                        if keys:
                            ### Keyword count                    
                            for prot in l:
                                jjj +=1
                                if jjj%1000==0:
                                    print "{:.1%}, {:7}/582865 ({:.2%} of the anotated) | cluster {:6}/{:6}".format(jjj/582865.,
                                                                                                                    jjj,
                                                                                                                    jj/lenanot,
                                                                                                                    C,
                                                                                                                    nb_clusters)
                                if pid(prot) in dataframe_prot.index:
                                    jj +=1
                                    func = str(g.data.proteins.ix[pid(prot),"function"])
                                    clusters["annotated"][-1] += 1
                                    for pos_key,query_key in queries:
                                        if func.find(query_key) != -1:
                                            key_matrix[pos_key,C] += 1
                                            key_count[pos_key] +=1
            if keys:
                with open(keyfile,"w") as f:
                    pickle.dump(key_matrix,f)
                    print("{} proteins, {} annotated.".format(jjj,jj,np.sum(clusters["annotated"])))
                    keywords["count"] = key_count 
                    store.append("keywords",keywords,format="table")
                    
            store.append("proteins",pandas.DataFrame(proteins).set_index("protein_id"), format="table")
            store.append("clusters",pandas.DataFrame(clusters).set_index("name"), format="table")
        else:
            try:
                with open(keyfile,"r") as f:
                    key_matrix = pickle.load(f)
            except Exception:
                key_matrix = None
        return store,key_matrix

