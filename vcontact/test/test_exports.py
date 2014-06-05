import numpy as np
import scipy.sparse as sparse
from .. import exports
import subprocess
import pandas 
import os 
import pandas.util.testing as ptest
import glob

F = {} # fixtures
def teardown():
    for f in glob.glob("testfile*"):
        print "removing {}".format(f)
        subprocess.call("head {}".format(f),shell=True)
        os.remove(f)

def setup():
    F["network"] = sparse.lil_matrix(np.matrix([[0,0,1,1,0,0],
                                                [0,0,0,0,0,1],
                                                [1,0,0,1,0,0],
                                                [1,1,1,0,0,0],
                                                [0,0,0,0,0,0], 
                                                [0,1,0,0,0,0]]))
    F["B"] = np.matrix([[1,0.1],
                        [0.1,1],
                        [1,0],
                        [1,0],
                        [0,0],
                        [0,1]],dtype=float)

    F["contigs"] =  pandas.DataFrame({"name":["c_{0}".format(i) for i in range(6)],
                                      "family": ["a","b","a","a","c","b"],
                                      "genus": ["a","b",None,"a","d","c"], # contig 2 has no known genus 
                                      "origin": ["test"]*6,
                                      "mcl_cluster": [0,1,0,0,np.nan,1],
                                      "membership":[1,1,1,1,0,1],
                                      "pos":range(6) })

def test_cytoscape_network():
    exports.cytoscape_network(F["network"],F["contigs"],"testfile",cluster=0,membership=F["B"])

