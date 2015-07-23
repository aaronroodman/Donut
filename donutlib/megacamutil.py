#
# $Rev:: 191                                                          $:  
# $Author:: roodman                                                   $:  
# $LastChangedDate:: 2014-09-03 11:00:33 -0700 (Wed, 03 Sep 2014)     $:
#
# Utility methods for MEGACam 4604x2048 pixels
# assuming 2x2 readout here for 2304x1024 pixels 13.5*2 = 27.0 microns
# plate scale = 169 micron/arcsec, or = 5.92 marcsec/micron 
# for 2x2 readout that is  = 0.160 arcsec/pixel
#
import numpy
from collections import OrderedDict

class megacaminfo(object):
    """ megacaminfo is a class used to contain Megacam geometry information and various utility routines
    """

    def info(self):
        """info returns a dictionary chock full of info on the Megacam geometry
        keyed by the CCD name
        """

        infoDict = OrderedDict()

        # store a dictionary for each CCD, keyed by the CCD name
        plateScale = 0.169 # mm/arcsec

        infoDict["IM1"] =  {"xCenter":  682.0*plateScale,"yCenter":-585*plateScale, "FAflag":False, "CCDNUM":1, "Offset": 0.0}
        infoDict["IM2"] =  {"xCenter":  511.0*plateScale,"yCenter":-585*plateScale, "FAflag":False, "CCDNUM":2, "Offset": 0.0}
        infoDict["IM3"] =  {"xCenter":  342.0*plateScale,"yCenter":-585*plateScale, "FAflag":False, "CCDNUM":3, "Offset": 0.0}
        infoDict["IM4"] =  {"xCenter":  171.0*plateScale,"yCenter":-585*plateScale, "FAflag":False, "CCDNUM":4, "Offset": 0.0}
        infoDict["IM5"] =  {"xCenter":    0.0*plateScale,"yCenter":-585*plateScale, "FAflag":False, "CCDNUM":5, "Offset": 0.0}
        infoDict["IM6"] =  {"xCenter": -171.0*plateScale,"yCenter":-585*plateScale, "FAflag":False, "CCDNUM":6, "Offset": 0.0}
        infoDict["IM7"] =  {"xCenter": -342.0*plateScale,"yCenter":-585*plateScale, "FAflag":False, "CCDNUM":7, "Offset": 0.0}
        infoDict["IM8"] =  {"xCenter": -511.0*plateScale,"yCenter":-585*plateScale, "FAflag":False, "CCDNUM":8, "Offset": 0.0}
        infoDict["IM9"] =  {"xCenter": -682.0*plateScale,"yCenter":-585*plateScale, "FAflag":False, "CCDNUM":9, "Offset": 0.0}

        infoDict["IM10"] =  {"xCenter":  682.0*plateScale,"yCenter":-185*plateScale, "FAflag":False, "CCDNUM":10, "Offset": 0.0}
        infoDict["IM11"] =  {"xCenter":  511.0*plateScale,"yCenter":-185*plateScale, "FAflag":False, "CCDNUM":11, "Offset": 0.0}
        infoDict["IM12"] =  {"xCenter":  342.0*plateScale,"yCenter":-185*plateScale, "FAflag":False, "CCDNUM":12, "Offset": 0.0}
        infoDict["IM13"] =  {"xCenter":  171.0*plateScale,"yCenter":-185*plateScale, "FAflag":False, "CCDNUM":13, "Offset": 0.0}
        infoDict["IM14"] =  {"xCenter":    0.0*plateScale,"yCenter":-185*plateScale, "FAflag":False, "CCDNUM":14, "Offset": 0.0}
        infoDict["IM15"] =  {"xCenter": -171.0*plateScale,"yCenter":-185*plateScale, "FAflag":False, "CCDNUM":15, "Offset": 0.0}
        infoDict["IM16"] =  {"xCenter": -342.0*plateScale,"yCenter":-185*plateScale, "FAflag":False, "CCDNUM":16, "Offset": 0.0}
        infoDict["IM17"] =  {"xCenter": -511.0*plateScale,"yCenter":-185*plateScale, "FAflag":False, "CCDNUM":17, "Offset": 0.0}
        infoDict["IM18"] =  {"xCenter": -682.0*plateScale,"yCenter":-185*plateScale, "FAflag":False, "CCDNUM":18, "Offset": 0.0}

        infoDict["IM19"] =  {"xCenter":  682.0*plateScale,"yCenter":185*plateScale, "FAflag":False, "CCDNUM":19, "Offset": 0.0}
        infoDict["IM20"] =  {"xCenter":  511.0*plateScale,"yCenter":185*plateScale, "FAflag":False, "CCDNUM":20, "Offset": 0.0}
        infoDict["IM21"] =  {"xCenter":  342.0*plateScale,"yCenter":185*plateScale, "FAflag":False, "CCDNUM":21, "Offset": 0.0}
        infoDict["IM22"] =  {"xCenter":  171.0*plateScale,"yCenter":185*plateScale, "FAflag":False, "CCDNUM":22, "Offset": 0.0}
        infoDict["IM23"] =  {"xCenter":    0.0*plateScale,"yCenter":185*plateScale, "FAflag":False, "CCDNUM":23, "Offset": 0.0}
        infoDict["IM24"] =  {"xCenter": -171.0*plateScale,"yCenter":185*plateScale, "FAflag":False, "CCDNUM":24, "Offset": 0.0}
        infoDict["IM25"] =  {"xCenter": -342.0*plateScale,"yCenter":185*plateScale, "FAflag":False, "CCDNUM":25, "Offset": 0.0}
        infoDict["IM26"] =  {"xCenter": -511.0*plateScale,"yCenter":185*plateScale, "FAflag":False, "CCDNUM":26, "Offset": 0.0}
        infoDict["IM27"] =  {"xCenter": -682.0*plateScale,"yCenter":185*plateScale, "FAflag":False, "CCDNUM":27, "Offset": 0.0}

        infoDict["IM28"] =  {"xCenter":  682.0*plateScale,"yCenter":585*plateScale, "FAflag":False, "CCDNUM":28, "Offset": 0.0}
        infoDict["IM29"] =  {"xCenter":  511.0*plateScale,"yCenter":585*plateScale, "FAflag":False, "CCDNUM":29, "Offset": 0.0}
        infoDict["IM30"] =  {"xCenter":  342.0*plateScale,"yCenter":585*plateScale, "FAflag":False, "CCDNUM":30, "Offset": 0.0}
        infoDict["IM31"] =  {"xCenter":  171.0*plateScale,"yCenter":585*plateScale, "FAflag":False, "CCDNUM":31, "Offset": 0.0}
        infoDict["IM32"] =  {"xCenter":    0.0*plateScale,"yCenter":585*plateScale, "FAflag":False, "CCDNUM":32, "Offset": 0.0}
        infoDict["IM33"] =  {"xCenter": -171.0*plateScale,"yCenter":585*plateScale, "FAflag":False, "CCDNUM":33, "Offset": 0.0}
        infoDict["IM34"] =  {"xCenter": -342.0*plateScale,"yCenter":585*plateScale, "FAflag":False, "CCDNUM":34, "Offset": 0.0}
        infoDict["IM35"] =  {"xCenter": -511.0*plateScale,"yCenter":585*plateScale, "FAflag":False, "CCDNUM":35, "Offset": 0.0}
        infoDict["IM36"] =  {"xCenter": -682.0*plateScale,"yCenter":585*plateScale, "FAflag":False, "CCDNUM":36, "Offset": 0.0}

        return infoDict
    

    def __init__(self,**inputDict):

        self.infoDict = self.info()

        nReadout = 2.0
        pixelSize = 0.0135 * nReadout
        self.mmperpixel = pixelSize
        self.rClear = 225.0

    def __getstate__(self):
        stateDict = {}
        keysToPickle = ['infoDict','mmperpixel','rClear']
        for key in keysToPickle:
            stateDict[key] = self.__dict__[key]
        return stateDict

    def __setstate__(self,state):
        for key in state:
            self.__dict__[key] = state[key]

    def getPosition(self,extname,ix,iy):
        """ return the x,y position in [mm] for a given CCD and pixel number
        note that the ix,iy are Image pixels - overscans removed - and start at zero
        """

        ccdinfo = self.infoDict[extname]

        # CCD size in pixels - assuming 2x2 operation
        xpixHalfSize = 512.
        ypixHalfSize = 1152.

        # calculate positions
        xPos = ccdinfo["xCenter"] + (float(ix)-xpixHalfSize+0.5)*self.mmperpixel
        yPos = ccdinfo["yCenter"] + (float(iy)-ypixHalfSize+0.5)*self.mmperpixel

        return xPos,yPos

    def getPixel(self,extname,xPos,yPos):
        """ given a coordinate in [mm], return pixel number
        """

        ccdinfo = self.infoDict[extname]

        # CCD size in pixels
        xpixHalfSize = 512.
        ypixHalfSize = 1152.

        # calculate positions
        ix = (xPos - ccdinfo["xCenter"]) / self.mmperpixel + xpixHalfSize - 0.5
        iy = (yPos - ccdinfo["yCenter"]) / self.mmperpixel + ypixHalfSize - 0.5

        return ix,iy

    def getSensor(self,xPos,yPos):
        """ given x,y position on the focal plane, return the sensor name
        or None, if not interior to a chip
        """
        for ext in self.infoDict.keys():
            ccdinfo = self.infoDict[ext]
            # is this x,y inside this chip?
            nxdif = numpy.abs( (xPos - ccdinfo["xCenter"]) / self.mmperpixel )
            nydif = numpy.abs( (yPos - ccdinfo["yCenter"]) / self.mmperpixel )

            # CCD size in pixels
            xpixHalfSize = 512.
            ypixHalfSize = 1152.

            if nxdif <= xpixHalfSize and nydif <= ypixHalfSize:
                return ext

        # get to here if we are not inside a chip
        return None
            
        

