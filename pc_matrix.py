""" Object containing th pc-profiles of the contigs and able to
    compute the similarity network on the contigs and the pcs. """
import logging
import sys
import pandas
import numpy as np
import scipy.stats as stats
import scipy.sparse as sparse
import networkx
import cPickle as pickle
logger = logging.getLogger(__name__)


class PCMatrix(object):
    """
    Protein Cluster presence/absence matrix, or PCprofiles.
       - Built from a protein cluster object.
       - Usefull in network building.
    """

    def __init__(self,cluster_proteins, clusters, ref_proteins, contigs, contig_names = None):
        """
        Args:
            cluster_proteins: pandas dataframe with :
                protein_id: (index)
                cluster:
            clusters: pandas dataframe with:
                name: (index)
                size:
            ref_proteins: = pandas dataframe with:
                protein_id: (index)
                contig:
                function:
            contigs: pandas dataframe with:
                name: (index)
                family:
                genus: 
                origin:

            contig_names: (iterable) names of contigs to keep in the profiles.
                (if none all contigs from contigs are used)
        """

        self.name = "PCmatrix"
        
        self.features = clusters
        self.features.index.name = "name"
        self.features.reset_index(inplace=True,drop=False)
        logger.debug(self.features.head())
        
        # Filtering the contigs 
        self.contigs = self.filtering(contigs, contig_names)
        logger.debug(self.contigs.head())

        # Building the PC-profiles
        self.matrix,self.singletons = self.load_profiles(ref_proteins, cluster_proteins)

        # Build the contig network.
        self.ntw = self.network(self.matrix,self.singletons)

        # Next step is probably to
        # Export a network to mcl : self.to_mcl(self.ntw)
        # Build the pc-network : self.pc_ntw = self.network_modules(self.matrix)


    def __repr__(self):
        return ("PC-profiles object {2} "
                "{0[0]} contigs by {0[1]} shared protein clusters "
                "and {1} singletons.").format(self.matrix.shape,
                                              self.singletons.sum(),
                                              self.name)


    def filtering(self,contigs, contig_names):
        """
        Filering the contigs 
        input :
        - contigs = pandas dataframe with:
              index : name, family, genus, origin        
        - contig_names = iterable : names of contigs to keep
          (if none all contigs from contigs are used)
        output 
        - contigs = filtered pandas dataframe 
        """
        # STEP 0 : Set index name.
        contigs.index.name = "name"
        contigs.reset_index(inplace=True,drop=False)

        # STEP 1 : filtering the contigs
        try:
            i = len(contigs)
            contigs = contigs.query("name in contig_names")
            logger.info("Filtered {} contigs".format(i-len(contigs)))
        except Exception as e:
            logger.debug("No filtering done (contig_name).")

        # Set the position of each contig in the matrix once and for all.
        contigs["pos"] = contigs.reset_index().index

        return contigs


    def load_profiles(self,ref_proteins,cluster_proteins):
        """
        Load a Protein cluster presence/absence matrix.

        INPUT:
        ref_proteins = pandas dataframe with:
              index : protein_id, contig, function

        cluster_proteins = pandas dataframe with :
              index : protein_id, cluster

        OUTPUT:
        - matrix (scipy.sparse.csr_matrix) contigs x protein_clusters.
          M(c,p) = Bool("Contig c has the PC p")
        """
        # A tale in three dataframe merging...

        ### 0 ### Associate each protein to its contig (position).
        # We perform an inner join to discard the proteins not belonging
        # our list of filtered contigs. 
        pc = pandas.merge(self.contigs, ref_proteins.reset_index(),
                          left_on="name", right_on="contig",
                          how="inner").loc[:,["pos","protein_id"]]
        pc.columns = ["pos_contig","protein_id"]
        #########
        #logger.debug(pc.head())

        ### 1 ### Associate each protein to its cluster (if any)
        # We join on the left so we keep all proteins. The ones without clusters have a NaN.
        pc = pandas.merge(pc,cluster_proteins.reset_index(),
                          left_on="protein_id",right_on="protein_id",
                          how="left").loc[:,["pos_contig", "cluster"]]
        ######### pc.columns are pos_contig (int), cluster (str)
        #logger.debug(pc.head())

        i = len(pc)
        singletons = pc.loc[pandas.isnull(pc.cluster)].groupby("pos_contig").pos_contig.count()
        pc.dropna(subset=["cluster"],inplace=True)
        #logger.debug("Droped {} singletons \n Prots :\n {} Singletons : \n {}".format(i-len(pc),pc.head(),singletons.head()))
        
        # Filter the protein cluster are not represented
        present_clusters = pc.cluster.dropna().drop_duplicates() 
        self.features = self.features.query("name in present_clusters")

        # Set the position of each feature in the matrix once and for all.
        self.features["pos"] = self.features.reset_index().index

        ### 2 ### Associate the clustrs to their position in the matrix.
        # Drop the duplicates as a contig need only to have one protein of the cluster.
        pc = pandas.merge(pc, self.features,
                          left_on="cluster", right_on="name",
                          how="left").loc[:,["pos_contig","pos"]].drop_duplicates()
        ######## pc.columns are pos_contig (int), pos (int, pos cluster)
        #logger.debug(pc.head())
        
        # Build the PC-profile matrix by putting ones in the coordinates given by pc
        matrix      = sparse.coo_matrix(([1]*len(pc),(zip(*pc.values))),dtype="bool")
        singletons = sparse.coo_matrix((singletons.values, (singletons.index.values,[0]*len(singletons)) ))
        matrix = sparse.csr_matrix(matrix)

        logger.info(("PC-profiles matrix {0[0]} contigs by {0[1]} shared protein clusters "
                      "and {1} singletons.").format(matrix.shape, singletons.sum()))
        return matrix,singletons

    def network(self,matrix,singletons,thres=1):
        """
        Compute the hypergeometric-similarity contig network.

        Args:
            matrix (scipy.sparse)x: contigs x protein clusters :
                M(c,p) == True <-> PC p is in Contig c.
            thres (float): Minimal significativity to store an edge value.
        
        Return
            scipy.sparse: S symmetric lil matrix, contigs x contigs.
          S(c,c) = sig(link)
        """

        # There are 
        contigs, pcs = matrix.shape
        pcs += singletons.sum()
        
        # Number of comparisons
        T = 0.5 * contigs * (contigs -1 )
        logT = np.log10(T)

        # Number of protein clusters in each contig
        # = # shared pcs + #singletons
        number_of_pc = matrix.sum(1) + singletons

        # Number of common protein clusters between two contigs
        commons_pc = matrix.dot(sparse.csr_matrix(matrix.transpose(),dtype=int))

        S = sparse.lil_matrix((contigs,contigs))

        # Display
        i = 0
        total_c = float(commons_pc.getnnz())

        for A,B in zip(*commons_pc.nonzero()) : #For A & B sharing contigs
            if A != B:
                # choose(a, k) * choose(C - a, b - k) / choose(C, b)
                # sf(k) = survival function = 1 -cdf(k) = 1 - P(x<k) = P(x>k)
                # sf(k-1)= P(x>k-1) = P(x>=k)
                # It is symmetric but I put the smallest before to avoid numerical bias.
                a,b = sorted([number_of_pc[A], number_of_pc[B]])
                pval = stats.hypergeom.sf(commons_pc[A,B]-1,pcs,a, b)
                sig = np.nan_to_num(-np.log10(pval)-logT)

                if sig>thres:
                    S[min(A,B),max(A,B)] = sig

                # Display
                i += 1
                if i%1000 == 0:
                    sys.stdout.write(".")
                if i%10000 == 0:
                    sys.stdout.write("{:6.2%} {}/{}\n".format(i/total_c,i,total_c))

        logger.debug("Hypergeometric similarity network : {0} genomes, {1} edges".format(contigs,S.getnnz()))
        S += S.T # Symmetry

        return S



    def network_modules(self,matrix=None,thres=1):
        """
        Compute the hypergeometric-similarity pc network.

        Warning:
            Use only the PCs that are present in 3 contigs or more.

        Args: 
            matrix (scipy.sparse): contigs x protein clusters :
                M(c,p) == True <-> PC p is in Contig c.
            thres (float): Minimal significativity to store an edge value.
        
        Returns:
            scipy.sparse: Symmetric lil_matrix, PCs x PCs
                S(c,c) = sig(link)
        """

        matrix = self.matrix if matrix == None else matrix
        contigs = matrix.shape[0]

        # Number of contig in which a given PC is found
        number_of_contigs = np.squeeze(np.asarray(np.transpose(matrix.sum(0))))

        # We only keep the pcs that are presents in more than 2 contigs
        pos_pcs_in_modules = [i for i,x in enumerate(number_of_contigs) if x>2]
        pcs_in_modules = len(pos_pcs_in_modules)
        logger.info("{} pc present in strictly more than 2 contigs".format(pcs_in_modules))

        # Filtering the matrix
        matrix = matrix[:,pos_pcs_in_modules]

        # Number of comparisons
        T = 0.5 * pcs_in_modules * (pcs_in_modules -1 )
        logT = np.log10(T)

        # Number of common contigs between two pcs
        commons_contigs = sparse.csr_matrix(matrix,dtype=int).transpose().dot(matrix)

        S = sparse.lil_matrix((pcs_in_modules,pcs_in_modules))
        i = 0
        total_c = float(commons_contigs.getnnz())
        for A,B in zip(*commons_contigs.nonzero()) :
            if A != B:
                # choose(a, k) * choose(C - a, b - k) / choose(C, b)
                # sf(k) = survival function = 1 -cdf(k) = 1 - P(x<k) = P(x>k)
                # sf(k-1)= P(x>k-1) = P(x>=k)
                # It is symmetric but I put the smallest before to avoid numerical biais.

                a,b = sorted([number_of_contigs[A], number_of_contigs[B]])
                pval = stats.hypergeom.sf(commons_contigs[A,B]-1,contigs,a, b)
                sig = np.nan_to_num(-np.log10(pval)-logT)


                if sig>thres:
                    S[min(A,B),max(A,B)] = sig

                i += 1
                if i%1000 == 0:
                    sys.stdout.write(".")
                if i%10000 == 0:
                    sys.stdout.write("{:6.2%} {}/{}\n".format(i/total_c,i,total_c))

        logger.debug("Hypergeometric similarity network : {0} pcs, {1} edges".format(pcs_in_modules,S.getnnz()))
        S += S.T # Symmetry

        return S

    def nodes_properties(self, matrix):
        """ Compute several node specific statistics.
        
        Args:
            matrix: (sparse.matrix) contig-similaritt network
        
        Side-effect:
            self.contigs: Add the following columns:
                degree: Number of edge to the node. 
                clustering_coefficient: Proportion of existing edges
                    over possible edges in the neighborhood of the node
                betweeness_centrality: sum of the fraction of all-pairs 
                    shortest path that pass trhough the node.
        """

        D = networkx.from_scipy_sparse_matrix(matrix)
        
        bc = pandas.Series(networkx.betweenness_centrality(D), name="betweeness_centrality")
        degr = pandas.Series(networkx.degree(D), name="degree")
        clcoef = pandas.Series(networkx.clustering(D), name="clustering_coef")
                
        df = pandas.concat([bc,degr,clcoef],axis=1)
        self.contigs = pandas.merge(self.contigs,df,left_on="pos",right_index=True)

    def to_pickle(self,path=None):
        """ Pickle (serialize) object to file path."""
        path = self.name+".pkle" if path is None else path
        with open(path, 'wb') as f:
            pickle.dump(self, f)

def read_pickle(path):
    """Read pickled object in file path."""
    with open(path, 'rb') as fh:
        return pickle.load(fh)


