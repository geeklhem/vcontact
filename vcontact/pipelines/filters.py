# coding: utf-8
import genome_clusters 
import pc_matrix 
import cPickle as pickle
import numpy as np
import pandas
import subprocess
import options
import exports


folder = options.data_folder + "tara/network/"
pelagiphages = ["NC_020481","NC_020482","NC_020483","NC_020484"]
names = ["HTVC010P", "HTVC011P", "HTVC019P", "HTVC008M"]
with open("{}tara/pc_matrix_object.pkle".format(options.data_folder),"r") as f:
    pcm = pickle.load(f)
thres = [1]#,5,10]

for t in thres :
    mcl_fi = "{}mcl20_sig{}_mcl20_contigs.clusters".format(folder,t)
    ntw_fi = "{}network_sig{}.ntw".format(folder,t)

    ntw_matrix = pcm.ntw.multiply(pcm.ntw>=t)
    pcm.to_mcl(ntw_matrix, ntw_fi)
    
    subprocess.call("mcl {0} --abc -I 2 -o {1}  ".format(ntw_fi, mcl_fi),shell=True)
    gc = genome_clusters.GenomeCluster((pcm.features,
                                        pcm.contigs.copy(),
                                        pcm.ntw.multiply(pcm.ntw>=t)),
                                       mcl_file=mcl_fi,
                                       name="tarac10mcl20_mcl20_sig{}".format(t))
    c_pel = [int(x) for x in gc.contigs.query("name in pelagiphages").mcl_cluster.drop_duplicates().values]

    
    gc.routine()
    gc.nodes_properties()
    gc.nodes_size()

    exports.cytoscape_network(gc.network,gc.contigs,
                              "{}network_sig_{}".format(folder,t),
                              cluster = c_pel ,
                              membership=gc.B)

    
    gc.contigs.to_pickle("{}contigs_sig{}.pandas".format(folder,t))
    gc.clusters.to_pickle("{}cluster_sig{}.pandas".format(folder,t))
    gc.summary.to_pickle("{}summary_sig{}.pandas".format(folder,t))
    

    with open("{}gc_object_sig{}.pkle".format(folder,t),"w") as f:
        pickle.dump(gc,f)


