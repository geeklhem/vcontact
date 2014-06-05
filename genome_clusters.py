"""Genome clusters : An object to work on the similarity network and
make the affiliations"""

import os
import subprocess
import cPickle as pickle
import logging

import numpy as np
import scipy.sparse as sparse
import pandas

import pc_matrix
import options
import matrices
import ml_functions
import associations

logger = logging.getLogger(__name__)

class GenomeCluster(object):
    """ Deal with the clusters of the contig similarity network

    Attributes:
            features: (pandas df) Protein clusters
            contigs: (pandas df) Contigs & Reference genomes
            network: (sparse matrix) Contig similarity network
            taxonomy: (pandas df) Taxonomic class
            clusters: (pandas df) Contig clusters
            mcl_results: (list of list) mcl_result[cluster][prot]
    """
    def __init__(self,pcm,inflation=2,threshold=None,name=None,membership_simple=False):
        """
        Init the object with a pc-profile object and perform the clustering

        Args:
            pcm: PCMatrix object or a tuple
                 (features (df), contigs (df), network (sparse matrix)).
            inflation (float): inflation for mcl.
            threshold (float): minimal significativity. 
            name (str): A name to identify the object.
            membership_simple (bool): if false use non boolean membership. 
        """
        self.name = "gc_sig{}_mcl{}".format(threshold,inflation) if name is None else name
        self.inflation = inflation
        self.thres = threshold

        if isinstance(pcm,pc_matrix.PCMatrix):
            self.features = pcm.features.copy()
            self.contigs = pcm.contigs.copy()
            self.network = pcm.ntw
        else:
            self.features,self.contigs,self.network = pcm[0].copy(), pcm[1].copy(), pcm.network

        if threshold is not None:
            before = self.network.getnnz()
            self.network = self.network.multiply(self.network>=self.thres)
            logger.debug(("Filtered {} edges according to the sig. threshold"
                           " {}.").format(before-self.network.getnnz(), self.thres))

            
        self.taxonomy = self.extract_taxonomy(self.contigs)

        self.clusters,self.mcl_results = self.clustering(options.data_folder+self.name,
                                                         self.inflation)


        self.matrix = {}
        logger.info("Computing membership matrix...")
        if membership_simple:
            self.matrix["B"] = matrices.bool_membership(self.contigs)
        else:
            self.matrix["B"] = matrices.membership(self.mcl_results,
                                                   self.network,
                                                   self.contigs)
        self.contigs = associations.contig_cluster(self.contigs,self.matrix["B"])


    def __repr__(self):
        return "GenomeCluster object {}, {} contigs and {} clusters".format(self.name,len(self.contigs),len(self.clusters))

    #--------------------------------------------------------------------------#
    # PART 1 : IMPORT, EXTRACTION, CLUSTERING
    #--------------------------------------------------------------------------#
    def extract_taxonomy(self, contigs, levels=("family","genus")):
        """ Build the taxonomy dataframe.

        Args:
            contigs (pandas.DataFrame): with colum "level" and "pos". 
            levels (list): column name in contigs

        Returns:
            dict: A dictionary of pandas.DataFrame, one key by taxonomic level.
        """
        contigs = self.contigs if contigs is None else contigs
        tax = {}
        for t in levels:
            tax[t] = pandas.DataFrame(contigs.groupby(t).pos.count(),
                                      columns=["references"])
            tax[t].index.name = "name"
            tax[t].reset_index(inplace=True)
            tax[t]["pos"] = tax[t].index
        return tax


    def clustering(self,basename, I, force=False):
        """Export the matrix, Run MCL and load the results

        Args:
            basename: (str) Path for the exported files
            I: (float) inflation for mcl
            force: (bool) overwrite existing file

        Returns:
            See self.load_clusters.

        Side-Effects:
           Save basename.ntw the network file
           Save basename.clusters the clustering results
           self.contig: add "cluster" column
        """

        fi_ntw = basename+".ntw"
        fi_clusters = basename+".clusters"

        # Export for MCL
        logger.info("Exporting for MCL")
        if not os.path.exists(fi_ntw) or force:
            self.to_mcl(self.network,fi_ntw)
        else:
            logger.debug("Network file already exist.")

        # MCL
        logger.info("Clustering the pc similarity-network")
        if not os.path.exists(fi_clusters or force):
            subprocess.call(("mcl {0} -o {1} --abc "
                             "-I {2}").format(fi_ntw,fi_clusters,I),
                            shell=True)
            logger.debug("MCL({}) results are saved in {}.".format(I,fi_clusters))
        else:
            logger.debug("MCL({}) file already exist.".format(I,fi_clusters))

        # Load clusters
        return self.load_clusters(fi_clusters)

    def to_mcl(self,matrix,fi,names=None):
        """Save a network in a file ready for MCL

        Args:
            matrix (scipy.sparse_matrix): network.
            fi (str): filename .
            names (pandas.dataframe): with the columns
                "pos":  (int) is the position in the matrix.
                "name": (str) column contain the name of the node.
                If None, self.contigs is used.

        Returns:
            str: filename
        """
        names = self.contigs if names == None else names
        names = names.set_index("pos").name
        with open(fi,"wb") as f:
            matrix = sparse.dok_matrix(matrix)
            for r,c in zip(*matrix.nonzero()):
                f.write(" ".join([str(x) for x in (names[r],
                                                   names[c],
                                                   matrix[r,c])]))
                f.write("\n")

        logger.debug("Saving network in file {0} ({1} lines).".format(fi,matrix.getnnz()))
        return fi

    def load_clusters(self,fi):
        """ Load clusters from the mcl results

        Args:
            fi (str): path to the MCL result file.

        Returns:
            df (pandas.DataFrame): give for each contig cluster
                its name, size and position in the matrix.

        Side-Effect:
            Modify self.contig to add the column "pos_cluster"
            giving the pos of the cluster it belongs to.

        The file fi was probably generated using :
        "mcl <file>.ntw --abc -I 2 -o <file>.clusters".
        """

        # Read the files
        with open(fi) as f:
            c = [ line.rstrip("\n").split("\t") for line in f ]
        c = [x for x in c if len(c)>1]
        name = ["cluster_{}".format(i) for i in range(len(c))]
        size = [len(i) for i in c]
        pos = range(len(c))

        logger.info(("{} clusters loaded (singletons and non-connected nodes "
                     "are dropped).").format(len(c)))

        # Update self.contigs (To refactor)
        self.contigs.reset_index(inplace=True)
        self.contigs.set_index("name",inplace=True)
        self.contigs["pos_cluster"] = np.nan

        for i,cluster in enumerate(c):
            for n in cluster:
                self.contigs.loc[n,"pos_cluster"] = i
        self.contigs.reset_index(inplace=True)

        return pandas.DataFrame({"name":name, "size":size,"pos":pos}),c


    #--------------------------------------------------------------------------#
    # PART 2: Affiliations
    #--------------------------------------------------------------------------#

    
    def total_affiliation(self,levels=("family","genus")):
        """Routine of the analysis using all the dataset.
        
        Args:
            levels (tuple): Taxonomic levels to consider.

        Returns:
            dataframe: Classification metrics. One line by taxonomic level.

        Warning:
            This function modify directly the attributes of the object.
        """

        results = []

        for level in levels:
            logger.info("Affiliation at the {} level...".format(level))
            self.matrix[level] = {}
            # Taxonomic reference matrix
            self.matrix[level]["K"] = matrices.reference_membership(level,
                                                                    self.contigs,
                                                                    self.taxonomy[level])

            # recall, precision and F-measure matrix
            (self.matrix[level]["Q"],
             self.matrix[level]["R"],
             self.matrix[level]["P"],
             self.matrix[level]["F"]) = matrices.correspondence(self.matrix[level]["K"],
                                                                self.matrix["B"])

            # Associate clusters and taxonomic classes.
            self.clusters, self.taxonomy[level] = associations.cluster_taxonomy(self.clusters,
                                                                                self.taxonomy[level],
                                                                                level,
                                                                                self.matrix[level]["P"],
                                                                                self.matrix[level]["R"])
            # Associate the contigs with the classes
            self.contigs =  associations.contig_taxonomy(self.contigs,
                                                         self.taxonomy[level],
                                                         self.clusters,
                                                         level)

            # Compute classification metrics.
            logger.info("Computing the classification metrics")
            results.append(ml_functions.classification_metrics(self.contigs.loc[:,[level,"predicted_"+level]],
                                                               ref_col=level,
                                                               pred_col="predicted_"+level))
        return pandas.DataFrame(results,levels)

    def cross_validation_affiliation(self,level="family", folds=10):
        """Cross validation affiliation.

        Cut the dataset (where a reference taxonomy exist) into <fold> equal parts
        and use alternativly <fold-2> of them to do the affiliation. Compute the
        classification metrics on the learning set and on the two remaining sets:
        cross validation (used to select the model) and test (used to have an 
        unbiased estimate of the error of the classification.)

        Note: 
            The splitting is stratified to keep the relative taxonomic classes in 
            the same proportions in each fold. 

        Args:
            level (str): Taxonomic level to consider.
            folds (int): number of folds for the cross-validation.

        Returns:
            dict: Dict of dataframes, one entry by set (learning,cv,test). 
                In the dataframes are the classification metrics, one row by
                selected cv-set. 
        """
        results = {"train_set":[],
                   "cv_set":[],
                   "test_set":[]}
        conditions = {}
        contigs = self.contigs
        taxonomy = self.taxonomy[level]
        clusters = self.clusters
        
        # Stratified split according to the taxonomic level. 
        contigs = ml_functions.split_dataset(contigs,level,folds)
        logger.info("{} folds cross-validation".format(folds))
        all_sets = frozenset(range(folds))

        for cv_set,test_set in zip(range(folds),range(1,folds)+[0]):
            logger.info(("Cross-validation fold {:2} "
                         "({:.0%})").format(cv_set,cv_set/float(folds)))

            # Filtering conditions: 
            train_set = list(all_sets - frozenset([cv_set,test_set]))
            conditions["cv_set"] = "cvset_{}=={}".format(level, cv_set)
            conditions["test_set"] = "cvset_{}=={}".format(level, test_set)
            conditions["train_set"] = "cvset_{} in {}".format(level, train_set)

            # Affiliation: 
            K = matrices.reference_membership(level, contigs,
                                              taxonomy, conditions["train_set"])
            _,R,P,_ =  matrices.correspondence(K, self.matrix["B"])
            clusters, taxonomy = associations.cluster_taxonomy(clusters, taxonomy,
                                                               level, P, R)
            contigs =  associations.contig_taxonomy(contigs, taxonomy, clusters,
                                                    level)

            # Computing metrics:
            for set_ in conditions.keys():
                df = contigs.query(conditions[set_]).loc[:,[level,"predicted_"+level]]
                results[set_].append(ml_functions.classification_metrics(df,
                                                                         ref_col=level,
                                                                         pred_col="predicted_"+level))

            # Cleaning the dataset for next step:
            contigs = contigs.drop("predicted_"+level,1)
            taxonomy = taxonomy.drop(["pos_cluster","recall"],1)
            clusters = clusters.drop(["pos_"+level,"precision_"+level],1)
        for set_ in results.keys():
            results[set_] = pandas.DataFrame(results[set_],[range(folds)])
        return results
    



    def learning_curve_affiliation(self,level="family", folds=10):
        """Learning curve affiliation.

        Cut the dataset (where a reference taxonomy exist) into <fold> equal parts
        and increasingly one to <fold-1> of them to do the affiliation. Compute the
        classification metrics on the learning set and on a randomly picked
        cross validation set in the non used subset. 

        Note: 
            The splitting is stratified to keep the relative taxonomic classes in 
            the same proportions in each fold. 

        Args:
            levels (str): Taxonomic level to consider.
            folds (int): number of folds for the cross-validation.

        Returns:
            dict: Dict of dict of dataframes, one entrie by level then one entry
                by set (learning,cv). In the dataframes are the classification
                metrics), one row by learning set size . 
        """
    
        results = {"train_set":[],
                   "cv_set":[]}
        conditions = {}
        contigs = self.contigs
        taxonomy = self.taxonomy[level]
        clusters = self.clusters
        
        # Stratified split according to the taxonomic level. 
        contigs = ml_functions.split_dataset(contigs,level,folds)
        logger.info("{} folds cross-validation".format(folds))
 

        for cv_set in range(folds):
            logger.info(("Cross-validation set {:2} "
                         "({:.0%})").format(cv_set,cv_set/float(folds)))

            remaining_sets = range(folds)
            remaining_sets.pop(cv_set)
            train_sets = [[remaining_sets[y] for y in range(x)] for x in range(1,folds)]
            for train_set in train_sets:
                logger.info("Training set of size {:2}".format(len(train_set)))

                # Filtering conditions: 
                conditions["cv_set"] = "cvset_{}=={}".format(level, cv_set)
                conditions["train_set"] = "cvset_{} in {}".format(level, train_set)

                # Affiliation: 
                K = matrices.reference_membership(level, contigs,
                                                  taxonomy, conditions["train_set"])
                _,R,P,_ =  matrices.correspondence(K, self.matrix["B"])
                clusters, taxonomy = associations.cluster_taxonomy(clusters, taxonomy,
                                                                   level, P, R)
                contigs =  associations.contig_taxonomy(contigs, taxonomy, clusters,
                                                        level)

                # Computing metrics:
                for set_ in results.keys():
                    df = contigs.query(conditions[set_]).loc[:,[level,"predicted_"+level]]
                    metrics = ml_functions.classification_metrics(df,
                                                                  ref_col=level,
                                                                  pred_col="predicted_"+level)
                    metrics["train_size"] = len(train_set)
                    metrics["cv_set"] = cv_set
                    results[set_].append(metrics)

                # Cleaning the dataset for next step:
                contigs = contigs.drop("predicted_"+level,1)
                taxonomy = taxonomy.drop(["pos_cluster","recall"],1)
                clusters = clusters.drop(["pos_"+level,"precision_"+level],1)
        for set_ in results.keys():
            results[set_] = pandas.DataFrame(results[set_])
        return results
    

    #--------------------------------------------------------------------------#
    # PART 4: Pickle-save
    #--------------------------------------------------------------------------#


    def to_pickle(self,path=None):
        """ Pickle (serialize) object to file path."""
        path = self.name+".pkle" if path is None else path
        with open(path, 'wb') as f:
            pickle.dump(self, f)

def read_pickle(path):
    """Read pickled object in file path."""
    with open(path, 'rb') as fh:
        return pickle.load(fh)
