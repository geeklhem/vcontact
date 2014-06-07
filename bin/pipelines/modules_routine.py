import logging
import cPickle as pickle
import os 
import scipy.sparse as sparse

#import numpy as np
import pandas
#import subprocess

import options
import modules
import pc_matrix 
#import genome_clusters
#import exports

#--------
pelagiphages = ["NC_020481","NC_020482","NC_020483","NC_020484"]
names = ["HTVC010P", "HTVC011P", "HTVC019P", "HTVC008M"]
#--------


logging.info("Loading pcm object")
with open("{}tara/pc_matrix_object.pkle".format(options.data_folder),"r") as f:
    pcm = pickle.load(f)    
folder = options.data_folder + "tara/network/"

if os.path.exists("{}/network_modules.pkle".format(folder)):
    logging.info("Loading the network of PCs")
    with open("{}/network_modules.pkle".format(folder),"r") as f:
        network_modules = pickle.load(f)    
else:
    logging.info("Building the network of PCs")
    network_modules = pcm.network_modules()
    with open("{}/network_modules.pkle".format(folder),"w") as f:
        pickle.dump(network_modules,f)    

logging.info("Loading contigs")
contigs = pandas.read_pickle(folder+"contigs_sig1.pandas")
clusters = pandas.read_pickle(folder+"cluster_sig1.pandas")

logging.info("Creating the Module object")
mod = modules.Modules((pcm.features,contigs,network_modules,pcm.matrix))

logging.info("Linking the modules with the clusters")
mod.matrix_modules_and_clusters = mod.link_modules_and_clusters(clusters)

with open("../data/cache/tara_contigs10_and_refseq_pblast_mcl20_keywords.pkle","r") as f:
   key_matrix = pickle.load(f)

c_pel = [int(x) for x in contigs.query("name in pelagiphages").mcl_cluster.drop_duplicates().values]
print modules.extract_modules(c_pel,mod,key_matrix)


pel_modules_data = modules.extract_modules(c_pel,mod,key_matrix)
pel_contigs = contigs.query("mcl_cluster in c_pel").sort(["mcl_cluster","origin","size"],
                                                         ascending=[True,False,False] )
pel_modules = pel_modules_data.pos_module

pel_pcs = mod.features.query("module in pel_modules").sort(["module","size"],ascending=[True,False])
pos_pel_pcs = pel_pcs.pos.values
pos_pel_contigs = pel_contigs.pos.values
pel_matrix = pcm.matrix[pos_pel_contigs,:][:,pos_pel_pcs] 
with open("../data/tara/modules/pelagiphages.pkle","w") as f:
    pickle.dump({"matrix":pel_matrix, "pcs":pel_pcs, "contigs":pel_contigs, "modules":pel_modules},f)