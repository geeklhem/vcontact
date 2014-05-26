"""Matrices used to compare the reference taxonomy with the clustering results.

B: Membership matrix.
    B[g,c]: proportion of edges weight linking node g to nodes in c.
K: Reference membership matrix.
    K[g,t]: (bool) node g belongs to taxonomic class t.
Q: Correspondance matrix.
    Q[c,t]: Sum for the sequence of class t of the membership to the clusters c. 
R: Recall matrix.
    R[c,t]: Proportion of (Q affiliated to t) that belongs to c.
P: Precision matrix.
    P[c,t]: Proportion of (Q in cluster c) that is affiliated to t.
F: F-measure matrix.
    F[c,t]: (2PR)/(P+R)
"""
import logging

import numpy as np
import scipy.sparse as sparse

logger = logging.getLogger(__name__)

def membership(mcl_results, network, nodes):
    """Membership matrix of the node to the clusters.

    Args:
        mcl_results: (list of list) a list of nodes name by cluster.
            Extracted from the lines of the mcl output.
        network: (sparse matrix) similarity network.
        nodes: (dataframe) with a column "name" corresponding to the
            entries in mcl_results and a column "pos" corresponding to
            the position in the matrix.

    Returns:
        sparse_matrix: B, Membership matrix #node X #clusters,
            B(g,c) is the proportion of edges weight linking the
            node g to the cluster C
    """

    network = network.todense()
    network = (network - network.min()) / (network.max() - network.min())

    #Sum of weights linking to a given target.
    B_sum = np.hstack([network.sum(1)]*len(mcl_results))


    # A list giving the columns of the i-th cluster members :
    clusters = [np.int_(nodes.query("name in members").ix[:, "pos"].values)
                for members in mcl_results]
    B_clust = np.hstack([network[:, cols].sum(1) for cols in clusters])


    B = B_clust/B_sum
    B = np.nan_to_num(B)
    return B


def reference_membership(level, contigs, taxonomy,
                         condition="origin=='refseq_jan14'"):
    """
    Build K the (Genome X Taxonomic class) reference membership matrix.

    Args:
        level: (str) column name of ref. taxonomy.
        contigs: (dataframe) with columns:
            pos: position in the matrix (int).
            "level": reference taxonomy (str).
        taxonomy: (dataframe) with columns:
            name: name of the class  (str).
            pos: position in matrix (int).
        condition: (str) a query to filter the contig dataframe.

    Returns:
        sparse_matrix: K, Bool(K[c,t]) == Contig c is of class t.
    """

    taxonomy = taxonomy.set_index("name")
    K = sparse.lil_matrix((len(contigs), len(taxonomy)), dtype=bool)
    contigs = contigs.query(condition).loc[:, ("pos", level)].dropna().values

    for pos, tax in contigs:
        K[pos, taxonomy.pos.loc[tax]] = 1

    K = sparse.csc_matrix(K)
    return K


def correspondence(K, B):
    """
    Build the (cluster X taxonomic class) correspondances matrix.

    Args:
        K (sparse matrix): Reference membership matrix.
        B (sparse matrix): Membership matrix.

    Returns:
        (tuple of sparse_matrix): including
            Q: Correspondance matrix
            R: Recall (or coverage)
            P: Precision
            F: F-measure 2PR/(P+R)
    """

    Q = np.dot(np.transpose(B), K.todense())

    # Precision
    Q_hsum = np.hstack([Q.sum(1)]*Q.shape[1])
    P = Q / Q_hsum
    P = np.nan_to_num(P)

    # Recall
    Q_vsum = np.vstack([Q.sum(0)]*Q.shape[0])
    R = Q / Q_vsum
    R = np.nan_to_num(R)

    F = 2 * np.divide(np.multiply(P, R), P+R)
    return Q, R, P, F

def clustering_wise_metrics(P, R, B, K):
    """Compute the clustering wise recall, precision and f-measure.

    Args:
        P (sparse matrix): Precision.
        R (sparse matrix): Recall.
        B (sparse matrix): Membership
        K (sparse matrix): Taxonomic classes

    Returns:
        tuple: including:
            int: Clustering wise precision
            int: Clustering wise recall
            int: Clustering wise f-measure
    """

    cwise_P = float(np.dot(B.sum(0), np.max(P, 1)))
    cwise_P /= B.sum()

    cwise_R = float(np.dot(np.max(R, 0), np.transpose(K.sum(0))))
    cwise_R /= K.sum()

    cwise_F = 2 * (cwise_P * cwise_R)/(cwise_P + cwise_R)
    return cwise_P, cwise_R, cwise_F
