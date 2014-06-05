import genome_clusters
import genomes 
import numpy as np
import pandas

def genus_by_clust(mcl_file="../data/tara/tara_c10_ntwk_mcl50.clusters"):
    """
    """
    print mcl_file
    contigs = genomes.Genomes(True,False,False).data.contigs
    with open(mcl_file) as f:
           clusters = [ line.rstrip("\n").split("\t") for line in f ]
        
    list_genus_by_clust = [contigs.query("name in members").ix[:,"genus"].drop_duplicates().dropna().values
                         for members in clusters]
    
    nb_clust = len(clusters)
    genus_by_clust_nb = [len(x) for x in list_genus_by_clust if len(x)]
    assoc = ["\t".join(x) for x in list_genus_by_clust if len(x)>1]
    index = [i for i,x in enumerate(list_genus_by_clust) if len(x)>1]
    for i in index:
        mm = clusters[i]
        print contigs.query("name in mm").dropna(subset=["genus"])
    size  = [len(c) for c,x in zip(clusters,list_genus_by_clust) if len(x)>1]
    print "\n".join(["{} (size {}): {}".format(x,s,y) for x,s,y in zip(index,size,assoc)]) 
    mean_genusbc = np.mean(genus_by_clust_nb)
    max_genusbc = np.max(genus_by_clust_nb)
    nb_genus_clust = len(list_genus_by_clust)
    aa = pandas.DataFrame(genus_by_clust_nb)
    aa.columns = ["nb"]
    print("Genus in {}/{} clusters - max genus by cluster {} (mean {})".format(nb_genus_clust,nb_clust,
                                                                               max_genusbc,mean_genusbc))
    print aa.groupby("nb").count()

    return genus_by_clust_nb


def family_by_clust(mcl_file="../data/tara/tara_c10_ntwk_mcl50.clusters"):
    """
    """
    print mcl_file
    contigs = genomes.Genomes(True,False,False).data.contigs
    with open(mcl_file) as f:
           clusters = [ line.rstrip("\n").split("\t") for line in f ]
        
    list_family_by_clust = [contigs.query("name in members").ix[:,"family"].drop_duplicates().dropna().values
                         for members in clusters]
    
    nb_clust = len(clusters)
    family_by_clust_nb = [len(x) for x in list_family_by_clust if len(x)]
    assoc = ["\t".join(x) for x in list_family_by_clust if len(x)>1]
    index = [i for i,x in enumerate(list_family_by_clust) if len(x)>1]
    for i in index:
        mm = clusters[i]
        print contigs.query("name in mm").dropna(subset=["family"])
    size  = [len(c) for c,x in zip(clusters,list_family_by_clust) if len(x)>1]
    print "\n".join(["{} (size {}): {}".format(x,s,y) for x,s,y in zip(index,size,assoc)]) 
    mean_familybc = np.mean(family_by_clust_nb)
    max_familybc = np.max(family_by_clust_nb)
    nb_family_clust = len(list_family_by_clust)
    aa = pandas.DataFrame(family_by_clust_nb)
    aa.columns = ["nb"]
    print("Family in {}/{} clusters - max family by cluster {} (mean {})".format(nb_family_clust,nb_clust,
                                                                               max_familybc,mean_familybc))
    print aa.groupby("nb").count()

    return family_by_clust_nb
