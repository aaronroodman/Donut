import numpy as np
import pickle
from scipy.interpolate import Rbf


class wavefrontmap(object):
    """ wavefrontmap is a class used to build and access a Wavefront map - zernike coefficients vs. X,Y

    Aaron Roodman (C) SLAC National Accelerator Laboratory, Stanford University 2018.
    """

    def __init__(self,file=""):
        # init contains all initializations which are done only once for all fits
        
        self.file = file
        mapdict = pickle.load(open(self.file,'rb'))
        self.x = mapdict['x']
        self.y = mapdict['y']
        self.zcoeff = mapdict['zcoeff']

        self.interpDict = {}
        for iZ in range(3,37):    # numbering is such that iZ=3 is zern4
            self.interpDict[iZ] = Rbf(self.x, self.y, self.zcoeff[:,iZ])        

    def get(self,x,y,nZernikeFirst=5,nZernikeLast=37):
        # fill an array with Zernike coefficients for this x,y in the Map

        zout = np.zeros((nZernikeLast-nZernikeFirst+1))
        for iZactual in range(nZernikeFirst,nZernikeLast+1):
            iZ = iZactual-1
            zout[iZactual-nZernikeFirst] = self.interpDict[iZ](x,y)

        return zout

    def getArray(self,x,y,iZactual):
        # x,y can be arrays and will return an array of Zernike values for coefficient iZactual
        iZ = iZactual - 1
        zout = self.interpDict[iZ](x,y)
        return zout
