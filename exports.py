from itertools import combinations
import numpy as np
import pandas 

def cytoscape_membership(B, contigs, fi):
    """Save the membership matrix in an input format for cytoscape (csv)
    INPUT:
    - genome_cluster (GenomesCluster) object
    - fi (str) filename (default self.name without parenthesis) 
    """

    groups = ["g{0}".format(n) for n in range(B.shape[1])]

    with open(fi,"wb") as f:
        f.write("Node\t")
        f.write("\t".join(groups))
        f.write("\tChart")
        f.write("\n")
        for r in range(B.shape[0]):
            f.write("{0}\t ".format(contigs.query("pos=={}".format(r)).name.values[0]))
            for c in range(B.shape[1]):
                f.write("{}\t".format(B[r,c]))
            f.write('piechart: attributelist="{0}" colorlist="contrasting"'.format(",".join(groups)))
            f.write("\n")


def cytoscape_network(network, dataframe, fi,
                      pos=None, cluster=None,
                      col_cluster="mcl_cluster", col_pos="pos", col_name="name",
                      membership=None):
    """Save the network defined by the matrix with
    the names given by dataframe into a file fi for cytoscape
    
    the position in the matrix and the node names are given by
    the pos_col and name_col resp.
    
    Save also all the informations about a node given by the dataframe
    in a info file"""

    fi_ntw = fi+"_cytoscape.ntw"
    fi_info = fi + "_cytoscape.info"
    fi_mb = fi + "_membership_cytoscape.info"
    fi_info_clust = fi + "_clusters_cytoscape.info"

    ### Get the position of the contigs in the selected clusters 
    dataframe = dataframe.reset_index().set_index(col_pos)
    other_clusters = list(dataframe.loc[:,col_cluster].drop_duplicates().dropna().values)
    if pos == None:
        if cluster != None:
                pos = [] 
                for c in cluster:
                        pos += list(dataframe.query("{}=={}".format(col_cluster,c)).index)
                        other_clusters.pop(other_clusters.index(c))
        else:
            raise ValueError("Need pos list or cluster index")
    print "{} contigs in those {} cluster(s), {} other clusters.".format(len(pos),len(cluster),len(other_clusters))

    ### Export the network 
    export_network(fi_ntw,network,dataframe,
                   pos,col_cluster,
                   other_clusters,membership) 

    ### Export the membership info 
    if membership != None :
        export_membership(fi_mb,membership,dataframe,pos)
    
    ### Export node informations
    dataframe = dataframe.loc[pos,:]
    dataframe.index.name = "pos2"
    dataframe = dataframe.reset_index().set_index(col_name)
    dataframe.to_csv(fi_info,sep="\t")
    print("Wrote {} lines in {}.".format(len(dataframe),fi_info))    

    ### Export clusters informations
    info_clust = pandas.DataFrame(dataframe.groupby("mcl_cluster").size(),columns=["size"])
    info_clust.to_csv(fi_info_clust,"\t")
    print("Wrote {} lines in {}.".format(len(info_clust),fi_info_clust))    

def export_membership(fi_mb, membership, dataframe, pos):
    """ Export the membership matrix """ 
    membership = membership[pos,:]
    groups = [n for n,i in enumerate((np.squeeze(np.asarray(membership.sum(0)))>0)) if i]
    membership = membership[:,groups]
    groups = ["g_{}".format(g) for g in groups]

    with open(fi_mb,"wb") as f:
        f.write("Node\t")
        f.write("\t".join(groups))
        f.write("\tChart")
        f.write("\n")
        line = 1
        for i,r in enumerate(pos):
                f.write("{0}\t ".format(dataframe.ix[r,"name"]))
                f.write("\t".join([str(i) for i in np.squeeze(np.asarray(membership[i,:]))]))

        f.write('\tpiechart: attributelist="{0}" showlabels=false colorlist="contrasting"'.format(",".join(groups)))
        f.write("\n")
        line += 1
    print("Wrote {} lines in {}.".format(membership.shape[0],fi_mb))    
    
    
def export_network(fi_ntw, network, dataframe,
                   pos, col_cluster,
                   other_clusters=None, membership=None):
    """ Export the network in a cytoscape file """ 
    line =1
    with open(fi_ntw,"w") as f:
        f.write("Node_1\t")
        f.write("Type\t")
        f.write("Node_2\t")
        f.write("Edge_weight")
        f.write("\n")
        for r,c in combinations(pos,2) :
                if network[r,c]:
                    n1 = dataframe.ix[r,"name"] 
                    n2 = dataframe.ix[c,"name"] 
                    f.write("\t".join([str(n1),
                                       "contig-contig",
                                       str(n2),
                                       str(network[r,c])]))
                    f.write("\n")
                    line +=1
        if other_clusters != None and membership != None:            
                for c in other_clusters:
                    pos_c = list(dataframe.query("{}=={}".format(col_cluster,c)).index)
                    if network[pos_c,:][:,pos_c].sum():
                        for n in pos:
                            s = network[n,pos_c].sum()
                            if s and membership[n,c] > 0.1:
                                f.write("\t".join([dataframe.ix[n,"name"],
                                                   "contig-cluster",
                                                   "cluster_{}".format(c),
                                                   str(s)]))
                                f.write("\n")
                                line += 1
    print("Wrote {} lines in {}.".format(line,fi_ntw))
