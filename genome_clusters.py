import pc_matrix
import numpy as np
import scipy.sparse as sparse
import logging
import pandas
import cPickle as pickle
class GenomeCluster(object):
    
    def __init__(self,pcm,mcl_file=None ):
        """
        pcm : PCMatrix object or tuble features (pandas df),contigs (pandas df),network (sparse matriex).
        """
        if isinstance(pcm,pc_matrix.PCMatrix):
            self.features = pcm.features
            self.contigs = pcm.contigs
            self.network = pcm.network
        else:
            self.features,self.contigs,self.network = pcm

        tax = []
        for t in ["family","genus"]:
            tax.append(pandas.DataFrame(self.contigs.loc[:,t]).drop_duplicates().dropna())
            tax[-1].columns = ["name"]
            tax[-1]["pos"] = range(len(tax[-1]))
            tax[-1]["level"] = t


        self.taxonomy = pandas.concat(tax)
        if mcl_file:
            self.clusters,self.mcl_results = self.load_clusters(mcl_file)


    def load_clusters(self,fi):
        """
        load clusters from the mcl results
        """
        with open(fi) as f:
           c = [ line.rstrip("\n").split("\t") for line in f ]
        name = ["cluster_{}".format(i) for i in range(len(c))]
        size = [len(i) for i in c]
        pos = range(len(c))
        
        return pandas.DataFrame({"name":name,"size":size,"pos":pos}),c 
        
    def membership_matrix(self,mcl_results,network=None,contigs=None,features=None):
        """Build the membership matrix
        INPUT:
        - mcl_results (list of list) a list of contigs name by cluster. Extracted from the lines of the mcl output
        - network (sparse matrix)
        OUTPUT:
        - (sparse matrix) Membership matrix #target X #clusters, B(g,c) is the proportion 
        of edges weight linking the node g to the cluster C"""

        # Parameters default values
        network = self.network if network == None else network
        contigs = self.contigs if contigs == None else contigs
        features = self.features if features == None else features
        
        network = (network - np.min(network)) / (np.max(network) - np.min(network))
        B_sum = np.hstack([network.sum(1)]*len(mcl_results)) #Sum of weights linking to a given target.

        
        # A list giving the columns of the i-th cluster members :
        logging.debug(contigs)
        clusters = [contigs.query("name in members").ix[:,"pos"].values for members in mcl_results]
        print clusters
        B_clust = np.hstack([network[:,cols].sum(1) for cols in clusters])


        B = B_clust/B_sum
        B = np.nan_to_num(B)
        return B

        
    def reference_membership_matrix(self,level,contigs=None):
        """
        Build K the (Genome X Taxonomic class) reference membership matrix. 
        INPUT :
        level (str) family or genus 
        """
        taxonomy = self.taxonomy.query("level=='{}'".format(level))
        contigs = self.contigs if contigs == None else contigs


        K = sparse.lil_matrix((len(contigs),len(taxonomy)), dtype=bool)
        for pos,tax in contigs.loc[:,("pos",level)].dropna().values:
            k = taxonomy.query("name==tax").pos.values 
            K[pos,k] = 1
        K = sparse.csc_matrix(K)
        return K

    def correspondence_matrix(self,K,B=None):
        """
        Build the (cluster X taxonomic class) correspondances matrix.  
        Arguments:
        - `K` (sparse matrix) Correspondance matrix
        - `B`: Membership matrix (default self.membership)  
        Returns:
        - Q : Correspondance matrix
        - R : Recall (or coverage)
        - P : Precision 
        """
        logging.info("Building Q R and P") 
                   
        
        Q = np.dot(np.transpose(B),  K.todense())

        # Precision
        Q_hsum = np.hstack([Q.sum(1)]*Q.shape[1])
        P = Q / Q_hsum
        P = np.nan_to_num(P)

        # Recall 
        Q_vsum = np.vstack([Q.sum(0)]*Q.shape[0]) 
        R = Q / Q_vsum
        R = np.nan_to_num(R)
                                    
        A = np.sqrt(np.array(P)*np.array(R))
        return Q,R,P,A

    def clustering_wise_pr(self,P,R,B,K):
        """Compute the clustering wise recall and precision
        """

        cwise_P = float(np.dot(B.sum(0),np.max(P,1)))
        cwise_P /= B.sum()

        cwise_R = float(np.dot(np.max(R,0),np.transpose(K.sum(0))))
        cwise_R /= K.sum()
        
        return cwise_P, cwise_R

    def associations(self,P,R,level):
        """
        Build associations
        Arguments:
        - `P`: Precision matrix (cluster X taxonomic Class) 
        - `R`: Recall matrix  (cluster X taxonomic Class) 
        - `lvl`: (str) taxonomic level of P & R 
        """

        df = {}
        df["pos"] = range(P.shape[0])
        df["pos_"+level] = np.squeeze(np.asarray(np.argmax(P,1)))
        df["precision_"+level] = np.squeeze(np.asarray(np.max(P,1)))
        df = pandas.DataFrame(df)
        df.loc[df["precision_"+level]==0,"pos_"+level] = np.nan

        self.clusters = pandas.merge(self.clusters,df,how='left')

        
        df = {}
        df["pos"] = range(R.shape[1])
        df["recall_"+level] = np.squeeze(np.asarray(np.max(R,0)))
        df["level"] = level
        df["pos_cluster_"+level] = np.squeeze(np.asarray(np.argmax(R,0)))
        df = pandas.DataFrame(df)
        df.loc[df["recall_"+level]==0,"pos_cluster_"+level] = np.nan

        self.taxonomy = pandas.merge(self.taxonomy,df,how='left')
     

    def affiliate(self,B=None):
        """
        Affiliate each contig with the maximal membership cluster.
        """
        cm = {}
        cm["pos_cluster"] = np.squeeze(np.asarray(np.argmax(B,1)))
        cm["membership"] = np.squeeze(np.asarray(np.max(B,1)))
        cm["pos"] = range(B.shape[0])
        cm = pandas.DataFrame(cm)
        cm.loc[cm["membership"]==0,"pos_cluster"] = np.nan
        self.contigs = pandas.merge(self.contigs,cm)

        
        aff =  pandas.merge(self.contigs,
                            self.clusters,
                            left_on="pos_cluster",
                            right_on="pos",
                            suffixes=["","_clusters"]).loc[:,( "name","family","genus",
                                                               "pos_family", "pos_genus", "membership")].sort(axis=1)
        for l in ("family","genus"):
            aff = pandas.merge(aff,self.taxonomy.query("level=='{}'".format(l)).loc[:,("pos","name")],
                               left_on="pos_{}".format(l), right_on="pos",
                               suffixes=["","_{}".format(l)])

        aff = aff.drop(["pos","pos_family","pos_genus"],1)
        aff.rename(columns={"name_family":"predicted_family",
                            "name_genus":"predicted_genus",
                            "family":"reference_family",
                            "genus":"reference_genus"},
                   inplace=True)

        return aff
    

