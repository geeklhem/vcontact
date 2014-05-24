"""
Modules are groups of protein families
"""

import pc_matrix
import numpy as np
import scipy.sparse as sparse
import logging
import pandas
import cPickle as pickle
import networkx
import options
from Bio import SeqIO
import subprocess
import os
import scipy.stats as stats
from itertools import combinations
logger = logging.getLogger(__name__)

class Modules(object):
    def __init__(self,pcm,threshold=10,name=None):
        """

        Args:
            pcm: PCMatrix object or a tuple:
                 (features (df), contigs (df),
                  pc similirity-network (sp.matrix), pc-profiles (sp.matrix)).
            threshold: (int) Minimal sig value to take into account an edge.
            name: (str) A name to identify the object.
        """

        self.thres = threshold
        self.name = "mod_sig{}".format(threshold) if name is None else name

        if isinstance(pcm,pc_matrix.PCMatrix):
            self.features = pcm.features
            self.contigs = pcm.contigs
            self.network = pcm.ntw_modules
            self.matrix = pcm.matrix
        else:
            self.features,self.contigs,self.network,self.matrix = pcm

        # Filter the network according to the threshold:
        before = self.network.getnnz()
        self.network = self.network.multiply(self.network>=self.thres)
        logger.debug(("Filtered {} edges according to the sig. threshold"
                       " {}.").format(before-self.network.getnnz(), self.thres))

        # Define the modules
        self.modules = self.define_modules(self.network)
        #self.matrix_module = self.module_in_contigs()
        #self.link_modules_and_clusters
    
    def __repr__(self):
        return ("Modules object {}, {} modules, (contigs, pc) : {},"
                "sig. threshold {} ").format(self.name, len(self.modules),
                                             self.matrix.shape, self.thres)

    def define_modules(self,matrix,I=5):
        """
        Save the pc network in a file ready for MCL
        Run MCL
        Load clusters from the mcl results

        Args:
            matrix: (scipy.sparse matrix) network.
            fi: (str) filename.
            I: (float) Inflation for mcl.

        Returns:
            A dataframe containing:
                name: the name of the module.
                size: the number of protein clusters in the module.
                pos: the position of the module in the matrix.
                proteins: the number of proteins in the module.
                annotated_proteins: the number of annotated proteins in the module.

        Side-Effects:
            self.features: Add the column "module".

        Saved Files:
            name.ntwk: The pc similarity network.
            name_mcl_I.clusters: mcl results.
            name_mcl_I_modules.pandas: the module dataframe.
            name_mcl_I_features.pandas: the pc dataframe.
        """

        basename = options.folders["modules"]+self.name
        fi_in = basename+".ntwk"
        fi_out = basename+"_mcl_{}.clusters".format(I)
        fi_dataframe = basename+"_mcl_{}_modules.pandas".format(I)
        fi_feat = basename+"_mcl_{}_features.pandas".format(I)


        # Save for MCL
        logger.info("Exporting the pc-network for mcl")
        names = self.features.set_index("pos").name 
        if not os.path.exists(fi_in):
            with open(fi_in,"wb") as f:
                matrix = sparse.dok_matrix(matrix)
                for r,c in zip(*matrix.nonzero()):
                    f.write(" ".join([str(x) for x in (names[r],
                                                       names[c],
                                                       matrix[r,c])]))
                    f.write("\n")

            logger.debug(("Saving network in file {0}"
                           " ({1} lines)").format(fi_in,matrix.getnnz()))
        else:
            logger.debug("Network file {} already exist.".format(fi_in))

        # Run MCL
        logger.info("Clustering the pc similarity-network")
        if not os.path.exists(fi_out):
            subprocess.call(("mcl {0} -o {1} --abc"
                             "-I {2}").format(fi_in,fi_out,I),
                            shell=True)
            logger.debug("MCL({}) results are saved in {}.".format(I,fi_out))
        else:
            logger.debug("MCL({}) file already exist.".format(I,fi_out))


        # Read the mcl results
        logger.info("Loading the clustering results")
        if not os.path.exists(fi_dataframe) or not os.path.exists(fi_feat):
            with open(fi_out) as f:
               c = [ line.rstrip("\n").split("\t")
                     for line in f]
            c = [x for x in c if len(x) > 1] # Drop singletons 
            name = ["module_{}".format(i) for i in range(len(c))]
            size = [len(i) for i in c]
            pos = range(len(c))
            proteins = np.zeros(len(c))
            annotated_proteins = np.zeros(len(c))

            self.features = self.features.reset_index().set_index("name")
            self.features["module"] = np.nan

            for i,cluster in enumerate(c):
                for n in cluster:
                    self.features.loc[n,"module"] = i
                    proteins[i] += self.features.loc[n,"size"]
                    annotated_proteins[i] += self.features.loc[n,"annotated"] 
            self.features.dropna(subset=["module"],inplace=True)
            
            
            dataframe = pandas.DataFrame({"name":name,
                                          "size":size,
                                          "pos":pos,
                                          "proteins":proteins,
                                          "annotated_proteins":annotated_proteins})
            dataframe.to_pickle(fi_dataframe)
            self.features.to_pickle(fi_feat)
            logger.debug(("Saving {} modules containing {} "
                          " protein clusters in {}.").format(len(name),
                                                             sum(size),
                                                             fi_dataframe))
        else:
            dataframe = pandas.read_pickle(fi_dataframe)
            self.features = pandas.read_pickle(fi_feat)
            logger.debug("Read {} modules from {}.".format(len(dataframe),
                                                            fi_dataframe))
        return dataframe


    def module_in_contigs(self,matrix=None,contigs=None,modules=None):
        """
        Compute the presence of modules in each contig

        Input :
            matrix: (scipy.sparse), M[contigs,pc] =  "pc in contig" (bool)
            contigs: (pandas.DataFrame)
            modules: (pandas.DataFrame) with columns : pos,

        Returns:
            S: (matrix), S[contig,module] = proportions of module's pcs
                present in contig (\in [0,1])
        """

        matrix = self.matrix if matrix is None else matrix
        contigs = self.contigs if contigs is None else contigs
        modules = self.modules if modules is None else modules

        matrix = matrix.tocsc()

        # Number of pcs of module m in each contig.
        N = sparse.lil_matrix((len(self.contigs),len(self.modules)))
        for m,data in self.features.reset_index().groupby("module"):
            pos = data.pos.values
            N[:,m] = matrix[:,pos].sum(1)

        # Number of pcs in each module
        sizes = np.matrix([self.modules.size.values]*N.shape[0])

        # Compute
        S = sparse.csc_matrix(N.todense()/sizes)

        return S

    def link_modules_and_clusters(self, clusters,
                                  matrix_modules=None,
                                  modules=None, contigs=None,
                                  thres=1,
                                  own_threshold=0.7):
        """
        Link the modules with the contigs clusters
        using the hypergeometric formula.

        Args:
            matrix_modules: (scipy.sparse)
                bool(M[contig,module]>=own_threshold): module is present in contig
            modules: (pandas.DataFrame)
            clusters: (pandas.DataFrame)
            thres: (float) significativity threshold.
            own_threshold: minimal proportion of PCs to "own" the module. 

        Returns:
            S (scipy.sparse) S[module,cluster] = sig.

        Formula:
            P(X>=c) ~ H(n,a,b)
                c: number of contigs in cluster c owning module m.
                n: number of contigs.
                a: number of contigs in the cluster c.
                b: number of contigs owning module m.
        """

        matrix_modules = self.matrix_module if matrix_modules is None else matrix_modules
        modules = self.modules if modules is None else modules
        contigs = self.contigs if contigs is None else contigs

        # Filtering 
        non_filtered = matrix_modules.getnnz()
        matrix_modules = matrix_modules >= own_threshold
        logger.info(("{} contigs-modules owning association, {} filtered "
                      "(a contig must have {:.0%} of the PCs to own a module)."
                      ).format(matrix_modules.getnnz(),
                               non_filtered-matrix_modules.getnnz(),
                               own_threshold))
        
        nb_contigs,nb_modules = matrix_modules.shape
        nb_clusters = len(clusters)
        matrix_module = matrix_modules.tocsc()
        logger.info(("Linking {} modules with "
                      "{} contigs clusters...").format(nb_modules,
                                                        nb_clusters))

        # Phage in a given cluster
        a_values = clusters.sort("pos").size.values

        # Phage displaying a given module
        b_values = np.squeeze(np.asarray(np.transpose(matrix_modules.sum(0))))


        # contig in cluster
        xy = contigs.reset_index().sort("pos").ix[:,["mcl_cluster","pos"]].dropna(subset=["mcl_cluster"]).values
        pa_matrix = sparse.coo_matrix( ([1]*len(xy), zip(*xy) ),
                                       shape=(nb_clusters,nb_contigs) )


        # Phage in a given cluster displaying a given module
        c_values = pa_matrix.dot(matrix_modules)

        # Number of comparisons
        logT = np.log10(nb_clusters*nb_modules)

        # Number of actual comparisons
        total_c = float(c_values.getnnz())

        S = sparse.lil_matrix((nb_clusters,nb_modules))
        for A,B in zip(*c_values.nonzero()) :
            # choose(a, k) * choose(C - a, b - k) / choose(C, b)
            # sf(k) = survival function = 1 -cdf(k) = 1 - P(x<k) = P(x>k)
            # sf(k-1)= P(x>k-1) = P(x>=k)
            pval = stats.hypergeom.sf(c_values[A,B]-1,nb_contigs,a_values[A], b_values[B])
            sig = np.nan_to_num(-np.log10(pval)-logT)
            if sig>thres:
                S[A,B] = sig
        logger.info("Network done {0[0]} clusters, {0[1]} modules and {1} edges.".format(S.shape,S.getnnz()))
        return S

    def link_modules_and_clusters_df(self,clusters,modules=None, matrix_modules_and_clusters=None):
        """Returns the dataframe giving the association between clusters and modules
        
        Args:
            clusters: (dataframe) informations about the clusters
            matrix_modules_and_clusters: (sparse matrix) significativity of the associations
            modules: (dataframe) information about the modules
        Returns:
            A dataframe containing the position of the module, 
            of the cluster, the significativity of the association 
            and merged with the lines from modules and clusters.. 
        """
        matrix = self.link_modules_and_clusters(clusters) if matrix_modules_and_clusters is None else matrix_modules_and_clusters
        modules = self.modules if modules is None else modules
        
        matrix = matrix.todok().items()
        data = pandas.DataFrame({'pos_module': [x[0][1] for x in matrix],
                                 "pos_cluster": [x[0][0] for x in matrix],
                                 "sig": [x[1] for x in matrix]})
        data = pandas.merge(data, modules,
                            left_on="pos_module",right_on="pos",
                            how="left")
        data = pandas.merge(data, clusters,
                            left_on="pos_cluster",right_on="pos",
                            how="left",
                            suffixes=["module","cluster"])
        return data

"""
def extract_modules(clusters,mod,keyword_matrix=None):
    matrix =  mod.matrix_modules_and_clusters[clusters,:].todok().items()
    pel_m_c = [x[0][0] for x in matrix]
    pel_m_s = [x[1] for x in matrix]
    pel_m = [x[0][1] for x in matrix]

    data = {'pos_module':pel_m,
            "pos_cluster":[],
            "sig":pel_m_s,
            "pcs":[],
            "proteins":[],
            "keys":[]}

    for s,m,c in zip(pel_m_s,pel_m,pel_m_c):
        data["pos_cluster"].append(clusters[c])
a = mod.features.query("module==m")
        data["pcs"].append(len(a))
        data["proteins"].append(a.size.sum())
        if keyword_matrix != None:
            count = keyword_matrix[:,a.pos.values].sum(1)
            data["keys"].append(", ".join(["{} ({})".format(k,int(count[n])) for n,k in enumerate(options.keywords) if count[n] ]))
        else:
            data["keys"].append(None)
    return pandas.DataFrame(data)
"""
