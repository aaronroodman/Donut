#
# $Rev:: 206                                                          $:  
# $Author:: roodman                                                   $:  
# $LastChangedDate:: 2015-08-03 11:19:38 -0700 (Mon, 03 Aug 2015)     $:
#
#
# Python Class to implement Inverse Distance Weighted interpolation
# using kNearestNeighbors.  
#
import numpy
from scipy.spatial import cKDTree
import pdb


class IDWInterp(object):
    """  IDWInterp implements an Inverse Distance Weighted interpolation.
    It may be also used for extrapolation, although the IDW method has no way of knowing
    how to extrapolate the distribution.  Results may be suspect near the Hull of
    the input point distributions, when the distribution varies quickly there.

    Usage:

    # X,Y,Z are 1D arrays with the input distribtions
    # kNN is the number of nearest neighbors to use
    # power is the weighting power, w =  1/Dist^power
    # epsilon is a added to all Distances - you should change the default!
    myIDW = IDWInterp(X,Y,Z,kNN=4,power=1)

    # returns Znew values at interpolant points atX,atY
    Znew = myIDW.ev(atX,atY)

    # IDW formula is
    # IDW formula is (for power=1)
    # value = (z1/d1 + z2/d2 + z3/d3) / (1/d1 + 1/d2 + 1/d3)

    """

    def __init__(self,X,Y,Z,kNN=4,power=1,epsilon=0.0):
        # initialize the class
        self.kNN = kNN
        self.power = power
        self.epsilon = epsilon
        self.X = X.copy()
        self.Y = Y.copy()
        self.Z = Z.copy()
        xyArr = numpy.array(zip(self.X,self.Y))

        # make the kNN tree
        self.kdtree= cKDTree(xyArr)

        # kNN to use
        if len(self.Z)<kNN:
            self.usekNN = len(self.Z)
            print "IDWInterp: asked for %d kNN, using %d" % (self.kNN,self.usekNN)
        else:
            self.usekNN = self.kNN

    def ev(self,X,Y):
        """ evaluate the interpolant
        """
        # find the nearest neighbors
        xyArr = numpy.array(zip(X,Y))
        dNeighbors,indexNeighbors = self.kdtree.query(xyArr,self.usekNN)

        #
        # get the Z values of the neighbors
        zNeighbors = self.Z[indexNeighbors]

        # add epsilon to all distances
        dNeighbors = dNeighbors + self.epsilon * numpy.ones(dNeighbors.shape)

        # IDW formula is (for power=1)
        # value = (z1/d1 + z2/d2 + z3/d3) / (1/d1 + 1/d2 + 1/d3)
        dwgtZArray = zNeighbors * (1.0/numpy.power(dNeighbors,self.power))
        dwgtArray = (1.0/numpy.power(dNeighbors,self.power))
        value = dwgtZArray.sum(1)/dwgtArray.sum(1)  # the (1) makes sure the sum is just over input entry

        # check them all for nan
#        for val in value:
#            if numpy.isnan(val):
#                print "IDWInterp: got a nan!?"
#                pdb.set_trace()
        return value
        
            
    
