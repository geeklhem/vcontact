""" Genome clusters : An object to work on the similarity network and make the affiliation"""
import pc_matrix
import numpy as np
import scipy.sparse as sparse
import logging
import pandas
import cPickle as pickle
import networkx
import options 
from Bio import SeqIO
logger = logging.getLogger(__name__)

class GenomeCluster(object):
    """ Genome cluster object """ 
    def __init__(self,pcm,mcl_file=None,name="gc"):
        """

        Args:
            pcm: PCMatrix object or a tuple
                 (features (df), contigs (df), network (sparse matrix)).
            mcl_file: (str) mcl results file.
            name: (str) A name to identify the object.
        """
        self.name = name
        
        if isinstance(pcm,pc_matrix.PCMatrix):
            self.features = pcm.features
            self.contigs = pcm.contigs
            self.network = pcm.ntw
        else:
            self.features,self.contigs,self.network = pcm

        self.taxonomy = self.extract_taxonomy(self.contigs)
        
        if mcl_file:
            self.clusters,self.mcl_results = self.load_clusters(mcl_file)

    def __repr__(self):
        return "GenomeCluster object {}, {} contigs and {} clusters".format(self.name,len(self.contigs),len(self.clusters))

    def extract_taxonomy(self, contigs, levels=["family","genus"]):
        """ Build the taxonomy dataframe.

        Args:
            contigs: (pandas.DataFrame) 
            levels: (list) 

        Returns:
            taxonomy: (pandas.DataFrame) 
        """
        tax = []
        for t in levels:
            tax.append(pandas.DataFrame(contigs.loc[:,t]).drop_duplicates().dropna())
            tax[-1].columns = ["name"]
            tax[-1]["pos"] = range(len(tax[-1]))
            tax[-1]["level"] = t
        taxonomy = pandas.concat(tax)
        return taxonomy

    def routine(self):
        """Routine of the analysis"""
        # Membership matrix
        self.B  = self.membership_matrix(self.mcl_results)

        # Taxonomic reference matrix
        self.Kf = self.reference_membership_matrix("family")
        self.Kg = self.reference_membership_matrix("genus")

        # recall, precision and accuracy matrix
        self.QRPA_f = self.correspondence_matrix(self.Kf,self.B)
        self.QRPA_g = self.correspondence_matrix(self.Kg,self.B)

        # Associate each cluster to a taxonomic class
        self.associations(self.QRPA_f[2],self.QRPA_f[1],"family")
        self.associations(self.QRPA_g[2],self.QRPA_g[1],"genus")

        # Associate each contig to its max mbship cluster
        self.aff = self.affiliate(self.B)

        # Summary 
        self.summary = self.compute_summary()
        
        logger.debug(self.summary)

    def load_clusters(self,fi):
        """ Load clusters from the mcl results

        Args: 
            fi: (str) MCL result file. 
        
        Returns:
            df: (pandas.DataFrame) give for each contig cluster 
                its name, size and position in the matrix.

        Side-Effect: 
            Modify self.contig to add the column "mcl_cluster"
            giving the pos of the cluster it belongs to. 

        The file fi was probably generated using :
        "mcl <file>.ntw --abc -I 2 -o <file>.clusters".
        """

        # Read the files
        with open(fi) as f:
           c = [ line.rstrip("\n").split("\t") for line in f ]
        name = ["cluster_{}".format(i) for i in range(len(c))]
        size = [len(i) for i in c]
        pos = range(len(c))

        # Update self.contigs (To refactor) 
        self.contigs.reset_index(inplace=True)
        self.contigs.set_index("name",inplace=True)
        self.contigs["mcl_cluster"] = np.nan

        for i,cluster in enumerate(c):
            for n in cluster:
                self.contigs.loc[n,"mcl_cluster"] = i
        self.contigs.reset_index(inplace=True)

        return pandas.DataFrame({"name":name, "size":size,"pos":pos}),c 
        
    def membership_matrix(self,mcl_results,network=None,contigs=None,features=None):
        """Build the membership matrix.

        Args:
            mcl_results: (list of list) a list of contigs name by cluster. 
                Extracted from the lines of the mcl output
            network: (sparse matrix)
        Returns:
            B: (sparse matrix) Membership matrix #target X #clusters, 
                B(g,c) is the proportion of edges weight linking the
                node g to the cluster C
        """

        # Parameters default values
        network = self.network   if network  is None else network
        contigs = self.contigs   if contigs  is None else contigs
        features = self.features if features is None else features

        network = network.todense()
        network = (network - network.min()) / (network.max() - network.min())

        B_sum = np.hstack([network.sum(1)]*len(mcl_results)) #Sum of weights linking to a given target.

        
        # A list giving the columns of the i-th cluster members :
        clusters = [np.int_(contigs.query("name in members").ix[:,"pos"].values) for members in mcl_results]
        B_clust = np.hstack([network[:,cols].sum(1) for cols in clusters])


        B = B_clust/B_sum
        B = np.nan_to_num(B)
        return B

    def link_clusters(self,mcl_results=None,network=None,contigs=None):
        """Link the clusters 
        
        Args:
            mcl_results: (list of list) a list of contigs name by cluster. 
                Extracted from the lines of the mcl output
            network: (sparse matrix)

        Returns:
            (sparse matrix) L matrix cluster*cluster, edges inter over edges intra. 
        """

        # Parameters default values
        mcl_results = self.mcl_results if mcl_results is None else mcl_results
        network = self.network if network is None else network
        contigs = self.contigs if contigs is None else contigs

        network = network.todense()
      
        network = (network - network.min()) / (network.max() - network.min())
      
        # A list giving the columns of the i-th cluster members :
        pos_nan = len(self.clusters)

        Z = sparse.coo_matrix(([1.0]*len(self.contigs),
                               [self.contigs.pos,self.contigs.pos_cluster.fillna(pos_nan)])).todense()
        


        S = np.transpose(Z).dot(network.dot(Z))
        

        sii = np.vstack( [S.diagonal()] *S.shape[0] )
        
        L = (S + S.T) / (sii + sii.T)
        
        np.fill_diagonal(L,0)
        L = np.nan_to_num(L)
        return L

    
    def reference_membership_matrix(self, level,
                                    contigs=None, taxonomy=None,
                                    reference_origin=["refseq_jan14","test"]):
        """
        Build K the (Genome X Taxonomic class) reference membership matrix. 
        
        Args:
            level: (str) family or genus.
            contigs: (dataframe) with columns:
                pos: position in the matrix (int).
                origin: used to select the contigs (str).
                family and genus: reference taxonomy (str). 
            taxonomy: (dataframe) with columns:
                level: family or genus (str).
                name: name of the class  (str).
                pos: position in matrix (int).
            reference_origin: (list) contigs to consider.

        Returns:
            K: (sparse matrix) Bool(K[c,t]) == Contig c is of class t.
        """
        taxonomy = self.taxonomy if taxonomy is None else taxonomy
        taxonomy = taxonomy.query("level=='{}'".format(level))
        contigs = self.contigs if contigs is None else contigs

        K = sparse.lil_matrix((len(contigs),len(taxonomy)), dtype=bool)
        for pos,tax,origin in contigs.loc[:,("pos",level,"origin")].dropna().values:
            if origin in reference_origin:
                k = taxonomy.query("name==tax").pos.values 
                K[pos,k] = 1
        K = sparse.csc_matrix(K)
        return K

    def correspondence_matrix(self,K,B=None):
        """
        Build the (cluster X taxonomic class) correspondances matrix.  

        Arguments:
            K: (sparse matrix) Correspondance matrix.
            B: (sparse matrix) Membership matrix (default self.membership).  

        Returns:
            Q: Correspondance matrix
            R: Recall (or coverage)
            P: Precision 
        """
        logger.info("Building Q R and P") 
        #print("Tried to multiply B : {} and K : {}".format(B.shape,K.shape))
        try:
            Q = np.dot(np.transpose(B),  K.todense())
        except ValueError as e:
            print("Tried to multiply B : {} and K : {}".format(B.shape,K.shape))
            raise
        
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

    def clustering_wise_pr(self, P, R, B, K):
        """Compute the clustering wise recall and precision.

        Args:
            P: (sparse matrix) Precision.
            R: (sparse matrix) Recall.
            B: (sparse matrix) Membership
            K: (sparse matrix) Taxonomic classes
        
        Returns:
            cwise_P: Clustering wise precision
            cwise_R: Clustering wise recall 
        """

        cwise_P = float(np.dot(B.sum(0),np.max(P,1)))
        cwise_P /= B.sum()

        cwise_R = float(np.dot(np.max(R,0),np.transpose(K.sum(0))))
        cwise_R /= K.sum()
        
        return cwise_P, cwise_R

    def associations(self, P, R, level):
        """Build associations

        Args:
            P: Precision matrix (cluster X taxonomic Class) 
            R: Recall matrix  (cluster X taxonomic Class) 
            lvl: (str) taxonomic level of P & R 

        Side-Effects:
            self.cluster: Add columns "pos_level" and "precision_level" 
            self.taxonomy: Add columns "pos_cluster_level" and "recall_level"
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
     

    def affiliate(self, B=None):
        """
        Affiliate each contig with the maximal membership cluster.

        Args:
            B: (scipy.sparse) Membership matrix
        
        Returns: 
            aff: (dataframe) 

        Side-effects: 
            self.contigs: 
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
                            suffixes=["","_clusters"]).loc[:,( "pos_cluster","name","family","genus","origin",
                                                               "pos_family", "pos_genus", "membership")].sort(axis=1)

        for l in ("family","genus"):
            aff = pandas.merge(aff,self.taxonomy.query("level=='{}'".format(l)).loc[:,("pos","name")],
                               how="left",
                               left_on="pos_{}".format(l), right_on="pos",
                               suffixes=["","_{}".format(l)])


        aff = aff.drop(["pos","pos_family","pos_genus"],1)
        aff.rename(columns={"name_family":"predicted_family",
                            "name_genus":"predicted_genus",
                            "family":"reference_family",
                            "genus":"reference_genus",
                            "pos_cluster":"cluster_max_membership"},
                   inplace=True)

        return aff

    def compute_summary(self):
        summary = {"clustering_wise_precision":[],
                   "clustering_wise_recall":[],
                   "level":[],
                   "name":[],
                   "classes":[],
                   "recall_micro":[],
                   "precision_micro":[],
                   "specificity_micro":[],
                   "accuracy_micro":[],
                   "recall_macro":[],
                   "precision_macro":[],
                   "specificity_macro":[],
                   "accuracy_macro":[],
                   "contigs":[],
                   "origin":[],
                   "affiliated_contigs":[],
                   "reference_contigs":[],
                   }

        for o in ["origin=='refseq_jan14'","origin!='refseq_jan14'"]:
            logger.info("FILTERING AFF for {}".format(o))
            filtered_aff = self.aff.query(o)
            if len(filtered_aff):
                for K,QRPA,level in zip([self.Kf,self.Kg], [self.QRPA_f,self.QRPA_g], ["family","genus"]):

                    summary["name"].append(self.name)
                    summary["level"].append(level)
                    summary["origin"].append(o)

                    cwp,cwr =  self.clustering_wise_pr(QRPA[2], QRPA[1], self.B, K)
                    summary["clustering_wise_precision"].append(cwp)
                    summary["clustering_wise_recall"].append(cwr)



                    classes = filtered_aff.loc[:,"reference_{}".format(level)].drop_duplicates().dropna().values

                    summary["classes"].append(len(classes))
                    summary["contigs"].append(len(self.contigs))


                    summary["recall_micro"].append(0)
                    summary["precision_micro"].append(0)
                    summary["specificity_micro"].append(0)
                    summary["accuracy_micro"].append(0)

                    conf = {"TP":0,"TN":0,"FP":0,"FN":0}
                    aff = filtered_aff.dropna(subset=["reference_{0}".format(level)])
                    summary["affiliated_contigs"].append(len(filtered_aff))
                    summary["reference_contigs"].append(len(aff)) 
                    for cat in classes :
                        TP = float(len(aff.query("reference_{0}=='{1}'&predicted_{0}=='{1}'".format(level,cat))))
                        TN = float(len(aff.query("reference_{0}!='{1}'&predicted_{0}!='{1}'".format(level,cat))))
                        FP = float(len(aff.query("reference_{0}!='{1}'&predicted_{0}=='{1}'".format(level,cat))))
                        FN = float(len(aff.query("reference_{0}=='{1}'&predicted_{0}!='{1}'".format(level,cat))))
                        conf["TP"] += TP
                        conf["TN"] += TN
                        conf["FP"] += FP
                        conf["FN"] += FN
                        if TP:
                            summary["recall_micro"][-1] += TP / (TP+FN)
                            summary["precision_micro"][-1] += TP / (TP+FP)
                        if TN:
                            summary["specificity_micro"][-1]+= TN / (FP + TN)
                        summary["accuracy_micro"][-1] += (TP + TN) / (TP+FP+FN+TN)
                        
                    if conf["TP"]+conf["FN"]:
                        summary["recall_macro"].append(conf["TP"]/(conf["TP"]+conf["FN"]))
                    else:
                        summary["recall_macro"].append(0)
                        
                    if conf["TP"]+conf["FP"]:
                        summary["precision_macro"].append(conf["TP"]/(conf["TP"]+conf["FP"]))
                    else:
                        summary["precision_macro"].append(0) 
                    if conf["FP"]+conf["TN"]:
                        summary["specificity_macro"].append(conf["TN"]/(conf["FP"]+conf["TN"]))
                    else:
                        summary["specificity_macro"].append(0)

                    if conf["TP"]+conf["TN"]+conf["FP"]+conf["FN"]:
                        summary["accuracy_macro"].append((conf["TP"]+conf["TN"])/(conf["TP"]+conf["TN"]+conf["FP"]+conf["FN"]))
                    else:
                        summary["accuracy_macro"].append(0)
                        
                    if len(classes):
                        summary["recall_micro"][-1] /= len(classes)
                        summary["precision_micro"][-1] /= len(classes)
                        summary["specificity_micro"][-1]/= len(classes)
                        summary["accuracy_micro"][-1] /= len(classes)
                    
                    
        return pandas.DataFrame(summary)

    def aff_contig(self,level="family",contigs=None,clusters=None,taxonomy=None):
        
        contigs = self.contigs if contigs is None else contigs
        clusters = self.clusters if clusters is None else clusters
        taxonomy = self.taxonomy if taxonomy is None else taxonomy
        
        m1 = pandas.merge(contigs.reset_index(),
                          clusters,
                          left_on="pos_cluster",right_on="pos",
                          how="left",
                          suffixes=["__contig","__cluster"])

        m2 = pandas.merge(m1,
                          taxonomy.query("level=='{}'".format(level)).loc[:,["pos","name"]],
                          left_on="pos_{}".format(level), right_on="pos",
                          how="left",
                          suffixes=["","__{}".format(level)])

        aff = m2.loc[:,["name__contig","name"]]
        aff.columns = ["contig",level]
        aff[level] = aff[level].fillna("Non affiliated")
        print aff.groupby(level).count()
        return aff  

    def nodes_properties(self):
        """ Compute several node specific statistics"""

        D = networkx.from_scipy_sparse_matrix(self.network)
        bc = networkx.betweenness_centrality(D)
        data_bc = pandas.DataFrame(pandas.Series(bc))
        data_bc.columns = ["betweeness_centrality"]
        self.contigs = pandas.merge(self.contigs,data_bc,left_on="pos",right_index=True)

        # Degree 
        degr = networkx.degree(D)
        degr = pandas.DataFrame(pandas.Series(degr))
        degr.columns = ["degree"]
        self.contigs = pandas.merge(self.contigs,degr,left_on="pos",right_index=True)

        # Clustering coeffciennt 
        clcoef = pandas.DataFrame(pandas.Series(networkx.clustering(D)))
        clcoef.columns = ["clustering_coef"]
        self.contigs = pandas.merge(self.contigs,clcoef,left_on="pos",right_index=True)

    def nodes_size(self):
        """Size in bp of the contigs. Should probably be moved in the genomes modules """
        if self.contigs.index.name != "name":
            self.contigs.reset_index(inplace=True)
            self.contigs.set_index("name",inplace=True)

  
        rec = SeqIO.parse(options.data_folder+"refseq/phage.genomic.gbff","gb")
        for r in rec:
            name = r.id.split(".")[0]
            self.contigs.loc[name,"size"] = len(r.seq)
            
        rec = SeqIO.parse(options.data_folder+"tara/tara_c10.fna","fasta")
        for r in rec:
            self.contigs.loc[r.id,"size"] = len(r.seq)

        self.contigs.sort("size",ascending=False,inplace=True)
        self.contigs["size_rank"] = [i for i,n in enumerate(self.contigs.size)]

    def cluster_network_cytoscape(self,fi_ntw=None,fi_clust_info=None,network=None):
        fi_ntw = self.name+"_contigsclusters_network.ntw"
        fi_clust_info = self.name + "_contigsclusters_network.info"
        info = self.clusters

        
        info = pandas.merge(info,self.taxonomy.query("level=='family'"),
                            how="left",
                            left_on="pos_family",right_on="pos",
                            suffixes=["","__family"]
                            )
        #print info
        info = pandas.merge(info,self.taxonomy.query("level=='genus'"),
                            how="left",
                            left_on="pos_genus",right_on="pos",
                            suffixes=["","__genus"]
                            )
        info.reset_index(inplace=True)
        info.set_index("name",inplace=True)
        #print info
        info = info.loc[:,["size","name__family","name__genus"]]
        info.to_csv(fi_clust_info,sep="\t")
       
        network = self.network if network is None else network
  

        L = self.link_clusters(network=network)
        
        cluster_idx = self.clusters.reset_index()
        cluster_idx.set_index("pos",inplace=True)
        #print cluster_idx
        with open(fi_ntw,"w") as f:
            f.write("Cluster_1\t")
            f.write("Cluster_2\t")
            f.write("inter_ov_intra")
            f.write("\n")
            for r in range(L.shape[0]):
                for c in range(r):
                    if L[r,c]:
                        n1 = cluster_idx.ix[r,"name"] if r != L.shape[0]-1 else "non-clustered" 
                        n2 = cluster_idx.ix[c,"name"] if c != L.shape[0]-1 else "non-clustered"
                        f.write("\t".join([str(n1),
                                           str(n2),
                                           str(L[r,c])]))
                        f.write("\n")        

    def to_pickle(self,path=None):
        """ Pickle (serialize) object to file path."""
        path = self.name+".pkle" if path is None else path  
        with open(path, 'wb') as f:
            pickle.dump(self, f)

def read_pickle(path):
    """Read pickled object in file path.""" 
    with open(path, 'rb') as fh:
        return pickle.load(fh)
