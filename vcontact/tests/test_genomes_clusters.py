import numpy as np
import scipy.sparse as sparse
from .. import genome_clusters
from .. import pc_matrix
from test_pc_matrix import setup
import pandas 
import os 
import pandas.util.testing as ptest

F = {} # fixtures
def teardown():
    os.remove("mcl_results")
    try:
        os.remove("gc_contigsclusters_network.ntw")
        os.remove("gc_contigsclusters_network.info")
    except OSError:
        print "File not found"
def setup():
    F["network_hypergeom"] = sparse.lil_matrix(np.matrix([[0,0,1,1,0,0],
                                                          [0,0,0,0,0,1],
                                                          [1,0,0,1,0,0],
                                                          [1,0,1,0,0,0],
                                                          [0,0,0,0,0,0], 
                                                          [0,1,0,0,0,0]]))
    F["network_hypergeom_interlink"] = sparse.lil_matrix(np.matrix([[0,1,1,1,0,0],
                                                                    [1,0,0,0,0,1],
                                                                    [1,0,0,1,0,0],
                                                                    [1,0,1,0,0,0],
                                                                    [0,0,0,0,0,0], 
                                                                    [0,1,0,0,0,0]]))

    
    F["B"] = np.matrix([[1,0],
                        [0,1],
                        [1,0],
                        [1,0],
                        [0,0],
                        [0,1]],dtype=float)

    with open("mcl_results","wb") as f:
        for l in [["c_0","c_2","c_3"], ["c_1","c_5"]]:
            f.write("\t".join(l)+"\n")

    F["contigs"] = pandas.DataFrame({"name":["c_{0}".format(i) for i in range(6)],
                                     "family": ["a","b","a","a","c","b"],
                                     "genus": ["a","b",None,"a","d","c"], # contig 2 has no known genus 
                                     "origin": ["test"]*6}).set_index("name")
    
    F["ref_proteins"] = pandas.DataFrame({"protein_id":["prot_{0}".format(i) for i in range(12)],
                                          "contig": ["c_0","c_0",
                                                     "c_1","c_1","c_1",
                                                     "c_2","c_2","c_2",
                                                     "c_3","c_3",
                                                     "c_5",
                                                     "c_0"], #prot not in the clustering
                                          "function": ["test"]*12}).set_index("protein_id")
     
    F["cluster_proteins"] =  pandas.DataFrame({"protein_id":["prot_{0}".format(i) for i in range(11)],
                                               "cluster": ["pc_{0}".format(i) for i in [1,2,
                                                                                       0,2,3,
                                                                                       0,1,2,
                                                                                       1,2,
                                                                                       3]]}).set_index("protein_id")
    
    F["clusters"] = pandas.DataFrame({"name":["pc_{0}".format(i) for i in range(4)],
                                      "size": [10,20,10,5]}).set_index("name")


    F["PCM"] = pc_matrix.PCMatrix(F["cluster_proteins"], F["clusters"],
                                  F["ref_proteins"], F["contigs"])

    F["gc"] = genome_clusters.GenomeCluster(F["PCM"],mcl_file="mcl_results")

    F["K_family"] = np.matrix([[1,0,0],
                               [0,1,0],
                               [1,0,0],
                               [1,0,0],
                               [0,0,1],
                               [0,1,0]])
    F["K_genus"] = np.matrix([[1,0,0,0],
                              [0,1,0,0],
                              [0,0,0,0],
                              [1,0,0,0],
                              [0,0,1,0],
                              [0,0,0,1]])
    F["Q_genus"] = np.matrix([[2,0,0,0],
                              [0,1,0,1]])
                                
    F["R_genus"] = np.matrix([[1,0,0,0],
                              [0,1,0,1]])
    
    F["P_genus"] = np.matrix([[1,0,0,0],
                              [0,0.5,0,0.5]])
    
    F["A_genus"] = np.matrix([[1,0,0,0],
                              [0,np.sqrt(0.5),0,np.sqrt(0.5)]])

    F["Q_family"] = np.matrix([[3,0,0],
                               [0,2,0]])
                                
    F["R_family"] = np.matrix([[1,0,0],
                               [0,1,0]])

    F["P_family"] = np.matrix([[1,0,0],
                               [0,1,0]])

    F["A_family"] = np.matrix([[1,0,0],
                               [0,1,0]])

    F["contigs_a"] = pandas.DataFrame({"name":["c_{0}".format(i) for i in range(6)],
                                       "family": ["a","b","a","a","c","b"],
                                       "genus": ["a","b",None,"a","d","c"], # contig 2 has no known genus 
                                       "origin": ["test"]*6,
                                       "pos_cluster": [0,1,0,0,np.nan,1],
                                       "membership":[1,1,1,1,0,1]})

    

    F["aff"] = pandas.DataFrame({"reference_family": ["a","b","a","a","b"],
                                 "reference_genus": ["a","b",None,"a","c"],
                                 "membership":[1,1,1,1,1],
                                 "origin":["test"]*5,
                                 "cluster_max_membership":[0,1,0,0,1],
                                 "name": ["c_0","c_1","c_2","c_3","c_5"],
                                 "predicted_family":["a","b","a","a","b",],
                                 "predicted_genus":["a","b","a","a","b",]})
    # --------| TP | TN | FP | FN | R | P | S | A 
    #---------|-----------------------------------
    # Fam   A |  3 | 2  | 0  | 0  | 1 | 1 | 1 | 1
    # Fam   B |  2 | 3  | 0  | 0  | 1 | 1 | 1 | 1
    #---------|-----------------------------------
    # Genus A |  2 | 2  | 0  | 0  | 1 | 1 | 1 | 1
    # Genus B |  1 | 2  | 1  | 0  | 1 | .5|2/3|3/4
    # Genus C |  0 | 3  | 0  | 1  | 0 | 0 | 1 |3/4
    # SUM ----|  3 | 7  | 1  | 1  |3/4|3/4|7/8|10/12
    
    F["summary"] = pandas.DataFrame({"clustering_wise_precision":[1,0.8],
                                     "clustering_wise_recall":[5/6.,0.8],
                                     "level":["family","genus"],
                                     "name":["gc","gc"],
                                     "classes":[2,3],
                                     "contigs":[6,6],
                                     "affiliated_contigs":[5,5],
                                     "reference_contigs":[5,4],
                                     "recall_micro":[1,2/3.0],
                                     "precision_micro":[1,1.5/3],
                                     "specificity_micro":[1,8/9.0],
                                     "accuracy_micro":[1,10/12.0],
                                     "recall_macro":[1,3/4.],
                                     "precision_macro":[1,3/4.],
                                     "specificity_macro":[1,7/8.],
                                     "accuracy_macro":[1,10/12.0],
                                     "origin":["origin!='refseq_jan14'"]*2,
                                 })
    

def test_membership_matrix():
    np.testing.assert_array_equal(F["gc"].membership_matrix(F["gc"].mcl_results,F["network_hypergeom"]),
                                      F["B"],
                                      "Membership matrix")
def test_ref_membership():
    for t in ("family","genus"):
        np.testing.assert_array_equal(F["gc"].reference_membership_matrix(t).todense(),
                                      F["K_"+t],
                                      "K ({} level)".format(t))

def test_correspondence_matrix():
    Q,R,P,A = F["gc"].correspondence_matrix(sparse.lil_matrix(F["K_genus"]), F["B"])
    np.testing.assert_array_equal(Q,F["Q_genus"],"Q_genus")
    np.testing.assert_array_equal(R,F["R_genus"],"R_genus")
    np.testing.assert_array_equal(P,F["P_genus"],"P_genus") 
    np.testing.assert_array_equal(A,F["A_genus"],"A_genus")

    Q,R,P,A = F["gc"].correspondence_matrix(sparse.lil_matrix(F["K_family"]), F["B"])
    np.testing.assert_array_equal(Q,F["Q_family"],"Q")
    np.testing.assert_array_equal(R,F["R_family"],"R")
    np.testing.assert_array_equal(P,F["P_family"],"P") 
    np.testing.assert_array_equal(A,F["A_family"],"A") 

def test_cwise_genus():
    pr = F["gc"].clustering_wise_pr(F["P_genus"],
                                    F["R_genus"],
                                    F["B"],
                                    sparse.lil_matrix(F["K_genus"]))
    np.testing.assert_array_almost_equal(pr,[0.8,0.8],2,"Genus cwise P and R")
    
    pr = F["gc"].clustering_wise_pr(F["P_family"],
                                    F["R_family"],
                                    F["B"],
                                    sparse.lil_matrix(F["K_family"]))
    np.testing.assert_array_almost_equal(pr,[1.0,0.833],2,"Family cwise P and R")

    

def test_affiliation():

    F["gc"].affiliate(F["B"])
    print F["gc"].contigs.sort(axis=1)
    print F["contigs_a"].sort(axis=1)
    ptest.assert_frame_equal(F["gc"].contigs.loc[:,('family', 'genus', 'membership',
                                                    'name', 'origin','pos_cluster',
                                                    )].sort(axis=1),F["contigs_a"].sort(axis=1),check_dtype=False)

def test_associations():
    F["gc"].associations(F["P_family"],F["R_family"],"family")
    F["gc"].associations(F["P_genus"],F["R_genus"],"genus")
    aff = F["gc"].affiliate(F["B"])
    print "\n\n"
    print "{:=^80}".format("OUTPUT")
    print aff.sort(axis=1).set_index("name").sort()
    print "\n\n"
    print "{:=^80}".format("REFERENCE")
    print F["aff"].sort(axis=1).set_index("name").sort()
    ptest.assert_frame_equal(aff.sort(axis=1).set_index("name").sort(),
                             F["aff"].sort(axis=1).set_index("name").sort(),
                             check_dtype=False)
    
def test_routine():
    F["gc"].routine()
    print "\n{:=^80}".format("REFERENCE")
    print F["summary"].sort(axis=1).set_index("level").sort(),

    print "\n{:=^80}".format("OUTPUT")
    print F["gc"].summary.sort(axis=1).set_index("level").sort()
    
    ptest.assert_frame_equal(F["gc"].summary.sort(axis=1).set_index("level").sort(),
                             F["summary"].sort(axis=1).set_index("level").sort(),
                             check_dtype=False)
    


def test_link():
    L =  F["gc"].link_clusters(network=F["network_hypergeom_interlink"])
    print L 

def test_cluster_cyto():
    F["gc"].cluster_network_cytoscape(network=F["network_hypergeom_interlink"])

    with open("gc_contigsclusters_network.ntw",'r') as f:
        print "\n".join([l for l in f])
    with open("gc_contigsclusters_network.info",'r') as f:
        print "\n".join([l for l in f])



