import numpy as np
import scipy.sparse as sparse
from .. import pc_matrix
import pandas 
import os 

F = {} # fixtures
def setup():
    a = np.array([[0,1,1,0],
                  [1,0,1,1],
                  [1,1,1,0],
                  [0,1,1,0],
                  [0,0,0,0],
                  [0,0,0,1]],dtype=bool)
    F["matrix"] = sparse.lil_matrix(a)
    F["network_hypergeom"] = np.matrix([[0,0,1,1,0,0],
                                        [0,0,0,0,0,1],
                                        [1,0,0,1,0,0],
                                        [1,0,1,0,0,0],
                                        [0,0,0,0,0,0], 
                                        [0,1,0,0,0,0]])
    F["contigs"] = pandas.DataFrame({"name":["c_{}".format(i) for i in range(6)],
                                     "family": ["a","b","a","a","c","b"],
                                     "genus": ["a","b",None,"a","d","c"], # contig 2 has no known genus 
                                     "origin": ["test"]*6}).set_index("name")
    F["ref_proteins"] = pandas.DataFrame({"protein_id":["prot_{}".format(i) for i in range(12)],
                                          "contig": ["c_0","c_0",
                                                     "c_1","c_1","c_1",
                                                     "c_2","c_2","c_2",
                                                     "c_3","c_3",
                                                     "c_5",
                                                     "c_0"], #prot not in the clustering
                                          "function": ["test"]*12}).set_index("protein_id")
    
    F["cluster_proteins"] =  pandas.DataFrame({"protein_id":["prot_{}".format(i) for i in range(11)],
                                               "cluster": ["pc_{}".format(i) for i in [1,2,
                                                                                       0,2,3,
                                                                                       0,1,2,
                                                                                       1,2,
                                                                                       3]]}).set_index("protein_id")
    
    F["clusters"] = pandas.DataFrame({"name":["pc_{}".format(i) for i in range(4)],
                                      "size": [10,20,10,5]}).set_index("name")

    F["PCM"] = pc_matrix.PCMatrix(F["cluster_proteins"], F["clusters"],
                                  F["ref_proteins"], F["contigs"])
                                    

def test_matrix():
    np.testing.assert_array_equal(F["PCM"].matrix.todense(), F["matrix"].todense())

def test_hypergeom():
    np.testing.assert_array_equal(F["PCM"].network(F["matrix"]).todense()>0, F["network_hypergeom"])
    F["PCM"].to_mcl( F["network_hypergeom"],"test.ntwk")

    
def teardown():
    os.remove("test.ntwk")

