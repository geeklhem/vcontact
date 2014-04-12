import pandas
import numpy as np
import scipy.stats as stats
import scipy.sparse as sparse
from itertools import combinations
import logging
import sys 
logging.basicConfig(level=logging.DEBUG)


class PCMatrix(object):
    """
    Protein Cluster presence/absence matrix, 
    - Built from a protein cluster object.
    - Usefull in network building, arff file creation.
    """
    
    def __init__(self,cluster_proteins,clusters,ref_proteins,contigs):
        """
        INPUT:
        cluster_proteins = pandas dataframe with :
              index : protein_id, cluster 
        
        clusters = pandas dataframe with:
              index : name, size

        ref_proteins = pandas dataframe with:
              index : protein_id, contig, function
        
        contigs = pandas dataframe with:
              index : name, family, genus, origin 
        """
        contigs.index.name = "name"
        self.contigs = contigs.reset_index()
        self.contigs["pos"] = self.contigs.index

        clusters.index.name = "name"
        self.features = clusters.reset_index() 
        self.features["pos"] = self.features.index

        #logging.debug("Contigs:\n {0}".format(self.contigs.head()))
        #logging.debug("Features:\n {0}".format(self.features.head()))
        #logging.debug("Proteins_clusters:\n {0}".format(cluster_proteins.head()))
        #logging.debug("Proteins_ref:\n {0}".format(ref_proteins.head()))

        self.matrix = self.load(ref_proteins,cluster_proteins)
        self.ntw = self.network(self.matrix)

    def load(self,ref_proteins,cluster_proteins):
        """Load a Protein cluster presence/absence matrix from a
        ProteinClusters object.
        INPUT :
        - df (pandas.Dataframe) with columns Genome and Cluster
        OUTPUT : 
        - matrix (sparse.matrix)
        - features (list)
        - contigs (pandas.Dataframe)
        """
        pc = ref_proteins.join(cluster_proteins,how="right")
        #logging.debug("Merge 0 :\n {0}".format(pc.head()))
        
        pc = pandas.merge(pc , self.features,
                          left_on="cluster",right_on="name").loc[:,["contig","pos"]].drop_duplicates()
    
        #logging.debug("Merge 1 :\n {0}".format(pc.head()))
        pc = pandas.merge(pc, self.contigs,
                          left_on="contig", right_on="name",
                          suffixes=["_cluster","_contig"] ).loc[:,["pos_contig","pos_cluster"]]

        #logging.debug("Merge 2 :\n {0}".format(pc.head()))

        matrix = sparse.coo_matrix(([1]*len(pc),(zip(*pc.values))),dtype="bool")
        matrix = sparse.csr_matrix(matrix)
        logging.debug("P/A Matrix {0[0]} contigs by {0[1]} protein clusters.".format(matrix.shape))
        return matrix 

    def network(self,matrix,thres=1):
        contigs, pcs = matrix.shape
        T = 0.5 * contigs * (contigs -1 )
        
        # Number of protein clusters in each contig 
        number_of_pc = matrix.sum(1)
        
        # Number of common protein clusters between two contigs 
        commons_pc = matrix.dot(sparse.csr_matrix(matrix.transpose(),dtype=int))
        

        
        S = sparse.lil_matrix((contigs,contigs))
        i = 0
        total_c = float(commons_pc.getnnz())
        for A,B in zip(*commons_pc.nonzero()) : # combinations(range(contigs),2):
            if A != B:
                # choose(a, k) * choose(C - a, b - k) / choose(C, b)
                # sf(k) = survival function = 1 -cdf(k) = 1 - P(x<k) = P(x>k) 
                pval = stats.hypergeom.sf(commons_pc[A,B],pcs, number_of_pc[A], number_of_pc[B]) 
                sig = np.nan_to_num(-np.log10(pval*T))

                if sig>thres:
                    S[min(A,B),max(A,B)] = sig
                i += 1
                if i%1000 == 0:
                    sys.stdout.write(".")
                if i%10000 == 0:  
                    sys.stdout.write("{:6.2%} {}/{}\n".format(i/total_c,i,total_c))
                    
        logging.debug("Hypergeometric similarity network : {0} genomes, {1} edges".format(contigs,S.getnnz()))
        S += S.T
        return S

    def to_arff(self,weka_file,class_value,matrix=None,feature_list=None,target=None):
        """ Convert a scipy sparse matrix (matrix) and a list of classess, to an arff weka file """
        pass
        """
        matrix = self.matrix if matrix == None else matrix
        feature_list = self.features if feature_list == None else feature_list
        target = self.target if target ==None else target 
        
        if matrix.shape == (1,1):
            logging.warning("File {0} not wrote. (Empty matrix).".format(weka_file))
            return 1

        if os.path.isfile(weka_file) and self.force:
            logging.info("File {0} not wrote. (File existing).".format(weka_file))
            return 0


        possible_classes = frozenset(class_value.values())
        relation_name = os.path.basename(weka_file).split(".")[0]
        logging.info("Relation name : {0} (file : {1})".format(relation_name,weka_file))

        if not sparse.isspmatrix_csr(matrix):
            matrix = sparse.csr_matrix(matrix)
        matrix.sort_indices()

        with open(weka_file,"wb") as f:
            f.write("@relation proteinclusters_{0}\n".format(relation_name))    
            for feature in feature_list:
                f.write("@attribute {0} numeric\n".format(feature))
            f.write('@attribute class {{"{0}"}}\n'.format('","'.join(possible_classes)))

            f.write("@data\n")
            for row_num in range(matrix.shape[0]):
                if target[row_num] in class_value :
                    cl = '"{0}"'.format(class_value[target[row_num]]) if target[row_num] in class_value else "?"
                    clusters = ["{0} {1}".format(c[0],c[1]) for c in zip(matrix[row_num,:].indices,matrix[row_num,:].data)]
                    f.write('{{{0}, {1} {2}}}\n'.format(", ".join(clusters),matrix.shape[1],cl))"""


    def to_mcl(self,matrix,fi):
        """Save a network in a file ready for MCL
        INPUT:
        - matrix (scipy.sparse matrix) network
        - fi (str) filename 
        """
               
        with open(fi,"wb") as f:
            matrix = sparse.dok_matrix(matrix)
            for r,c in zip(*matrix.nonzero()):
                f.write(" ".join([str(x) for x in (self.contigs.ix[r,"name"],
                                                   self.contigs.ix[c,"name"],
                                                   matrix[r,c])]))
                f.write("\n")

        logging.debug("Saving network in file {0} ({1} lines).".format(fi,matrix.getnnz()))
        return fi


    def to_cytoscape():
        """
        """
        pass 
