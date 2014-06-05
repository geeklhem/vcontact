""" test feature select"""

from .. import feature_select
import numpy as np 
F = {}
def setup():
    F["vectors"] = np.array([[3,2],
                              [0,1]])
    F["gram_smithed"] = np.array([[3,0],
                                   [0,1]])


def test_gram_smithed():
    print "\nINPUT : \n{}".format(F["vectors"])
    print "\nEXPECTED : \n {}".format(F["gram_smithed"])
    
    out = feature_select.gram_smith(F["vectors"])
    
    print "\n OUTPUT \n {}".format(out)

    np.testing.assert_array_equal(out,F["gram_smithed"])

