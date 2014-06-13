""" Object containing the pc-profiles of the contigs and able to
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


class PCProfiles(object):
    """
    Protein Cluster presence/absence matrix, or PCprofiles.
    Used to build the similarity networks.

    Attributes:
        contigs (pandas.DataFrame):
        pcs (pandas.DataFrame):
        name (str): 
        profiles (sparse.matrix):
        singletons (sparse.matrix): 
        contig_ntw (sparse.matrix):
        modules_ntw (sparse.matrix):
    """

    def __init__(self,contigs,pcs,profiles,name=None,
                 sig=1.0,sig_mod=1.0,mod_shared_min=3):
        """
        Args:
            contigs (dataframe): Contig info, required fields are pos, id.
            pcs (dataframe): Protein clusters info, required fields are pos, id.
            profiles (dict): Required field are matrix (the contigXpc profiles
                matrix) and singletons (the contigsX1 matrix giving the number
                of singletons by contig.)
            sig (float): Sig. threshold in the contig similarity network.
            sig_mod (float): Sig. threshold in the pc similarity network.
            mod_shared_min (float): Minimal number of contigs a pc must appeear into
                to be taken into account in the modules computing.
            name (str): name the object (useful in interactive mode)
        """

        
        if name is None:
            self.name = "PCprofiles"

        # Get the data
        self.contigs = contigs 
        self.pcs = pcs 
        self.matrix,self.singletons = profiles["matrix"], profiles["singletons"]

        # Store the parameters 
        self.sig = sig
        self.sig_mod = sig_mod
        self.mod_shared_min = mod_shared_min
        
        
        # Copute the networks 
        self.ntw = self.network(self.matrix,self.singletons,thres=sig)
        self.ntw_modules = self.network_modules(self.matrix,
                                                thres=sig_mod,
                                                mod_shared_min=mod_shared_min)

    def __repr__(self):
        return ("PC-profiles object {2} "
                "{0[0]} contigs by {0[1]} shared protein clusters "
                "and {1} singletons.\n"
                "Contig sig threshold: {}, Module sig threshold: {}"
                "Module shared min {}").format(self.matrix.shape,
                                               self.singletons.sum(),
                                               self.name, self.sig,
                                               self.sig_mod, self.mod_shared_min)


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
        number_of_pc = number_of_pc.A1 #Transform into a flat array

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
                sig = min(300,np.nan_to_num(-np.log10(pval)-logT))
                if sig>thres:
                    S[min(A,B),max(A,B)] = sig
                # Display
                i += 1
                if i%1000 == 0:
                    sys.stdout.write(".")
                if i%10000 == 0:
                    sys.stdout.write("{:6.2%} {}/{}\n".format(i/total_c,i,total_c))
        S += S.T # Symmetry
        S = S.tocsr()
        if len(S.data) != 0:
            logger.debug(("Hypergeometric contig-similarity network : {0} contigs, "
                          "{1} edges (min:{2:.2} max: {3:.2}, threshold was {4})").format(contigs,
                                                                                          S.getnnz(),
                                                                                          S.data.min(),
                                                                                          S.data.max(),
                                                                                          thres))
        else:
            raise ValueError("No edge in the similarity network !") 
        return S



    def network_modules(self,matrix=None,thres=1,mod_shared_min=3):
        """
        Compute the hypergeometric-similarity pc network.

        Warning:
            Use only the PCs that are present in 3 contigs or more.

        Args: 
            matrix (scipy.sparse): contigs x protein clusters :
                M(c,p) == True <-> PC p is in Contig c.
            thres (float): Minimal significativity to store an edge value.
            mod_shared_min (float): Minimal number of contigs a pc must appeear into
                to be taken into account in the modules computing.

        
        Returns:
            scipy.sparse: Symmetric lil_matrix, PCs x PCs
                S(c,c) = sig(link)
        """

        matrix = self.matrix if matrix == None else matrix
        contigs = matrix.shape[0]

        # Number of contig in which a given PC is found
        number_of_contigs = np.squeeze(np.asarray(np.transpose(matrix.sum(0))))

        # We only keep the pcs that are presents in more than "mod_shared_min" contigs
        pos_pcs_in_modules = [i for i,x in enumerate(number_of_contigs) if x>=mod_shared_min]
        pcs_in_modules = len(pos_pcs_in_modules)
        logger.debug("{} pc present in strictly more than {} contigs".format(pcs_in_modules,mod_shared_min))

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
                sig = min(300,np.nan_to_num(-np.log10(pval)-logT))


                if sig>thres:
                    S[min(A,B),max(A,B)] = sig

                i += 1
                if i%1000 == 0:
                    sys.stdout.write(".")
                if i%10000 == 0:
                    sys.stdout.write("{:6.2%} {}/{}\n".format(i/total_c,i,total_c))

        logger.debug("Hypergeometric pcs-similarity network : {0} pcs, {1} edges".format(pcs_in_modules,S.getnnz()))
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


