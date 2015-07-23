#
# $Rev:: 191                                                          $:  
# $Author:: roodman                                                   $:  
# $LastChangedDate:: 2014-09-03 11:00:33 -0700 (Wed, 03 Sep 2014)     $:
#
# Utility methods for DECam
#
import numpy
from collections import OrderedDict

class decaminfo(object):
    """ decaminfo is a class used to contain DECam geometry information and various utility routines
    """

    def info(self):
        """info returns a dictionary chock full of info on the DECam geometry
        keyed by the CCD name
        """

        infoDict = OrderedDict()

        # store a dictionary for each CCD, keyed by the CCD name
        # AJR 9/14/2012 fixed these to agree with the DS9 coordinate system
        infoDict["S1"] =  {"xCenter":  -16.908,"yCenter":-191.670, "FAflag":False, "CCDNUM":25, "Offset": 0.0}
        infoDict["S2"]  = {"xCenter":  -16.908,"yCenter":-127.780, "FAflag":False, "CCDNUM":26, "Offset": 0.0}
        infoDict["S3"]  = {"xCenter":  -16.908,"yCenter": -63.890, "FAflag":False, "CCDNUM":27, "Offset": 0.0} 
        infoDict["S4"]  = {"xCenter":  -16.908,"yCenter":   0.000, "FAflag":False, "CCDNUM":28, "Offset": 0.0}
        infoDict["S5"]  = {"xCenter":  -16.908,"yCenter":  63.890, "FAflag":False, "CCDNUM":29, "Offset": 0.0}
        infoDict["S6"]  = {"xCenter":  -16.908,"yCenter": 127.780, "FAflag":False, "CCDNUM":30, "Offset": 0.0}
        infoDict["S7"]  = {"xCenter":  -16.908,"yCenter": 191.670, "FAflag":False, "CCDNUM":31, "Offset": 0.0}
        infoDict["S8"]  = {"xCenter":  -50.724,"yCenter":-159.725, "FAflag":False, "CCDNUM":19, "Offset": 0.0}
        infoDict["S9"]  = {"xCenter":  -50.724,"yCenter": -95.835, "FAflag":False, "CCDNUM":20, "Offset": 0.0}
        infoDict["S10"] = {"xCenter":  -50.724,"yCenter": -31.945, "FAflag":False, "CCDNUM":21, "Offset": 0.0}
        infoDict["S11"] = {"xCenter":  -50.724,"yCenter":  31.945, "FAflag":False, "CCDNUM":22, "Offset": 0.0}
        infoDict["S12"] = {"xCenter":  -50.724,"yCenter":  95.835, "FAflag":False, "CCDNUM":23, "Offset": 0.0}
        infoDict["S13"] = {"xCenter":  -50.724,"yCenter": 159.725, "FAflag":False, "CCDNUM":24, "Offset": 0.0}
        infoDict["S14"] = {"xCenter":  -84.540,"yCenter":-159.725, "FAflag":False, "CCDNUM":13, "Offset": 0.0}
        infoDict["S15"] = {"xCenter":  -84.540,"yCenter": -95.835, "FAflag":False, "CCDNUM":14, "Offset": 0.0}
        infoDict["S16"] = {"xCenter":  -84.540,"yCenter": -31.945, "FAflag":False, "CCDNUM":15, "Offset": 0.0}
        infoDict["S17"] = {"xCenter":  -84.540,"yCenter":  31.945, "FAflag":False, "CCDNUM":16, "Offset": 0.0}
        infoDict["S18"] = {"xCenter":  -84.540,"yCenter":  95.835, "FAflag":False, "CCDNUM":17, "Offset": 0.0}
        infoDict["S19"] = {"xCenter":  -84.540,"yCenter": 159.725, "FAflag":False, "CCDNUM":18, "Offset": 0.0}
        infoDict["S20"] = {"xCenter": -118.356,"yCenter":-127.780, "FAflag":False, "CCDNUM":8 , "Offset": 0.0}
        infoDict["S21"] = {"xCenter": -118.356,"yCenter": -63.890, "FAflag":False, "CCDNUM":9 , "Offset": 0.0}
        infoDict["S22"] = {"xCenter": -118.356,"yCenter":   0.000, "FAflag":False, "CCDNUM":10, "Offset": 0.0}
        infoDict["S23"] = {"xCenter": -118.356,"yCenter":  63.890, "FAflag":False, "CCDNUM":11, "Offset": 0.0}
        infoDict["S24"] = {"xCenter": -118.356,"yCenter": 127.780, "FAflag":False, "CCDNUM":12, "Offset": 0.0}
        infoDict["S25"] = {"xCenter": -152.172,"yCenter": -95.835, "FAflag":False, "CCDNUM":4 , "Offset": 0.0}
        infoDict["S26"] = {"xCenter": -152.172,"yCenter": -31.945, "FAflag":False, "CCDNUM":5 , "Offset": 0.0}
        infoDict["S27"] = {"xCenter": -152.172,"yCenter":  31.945, "FAflag":False, "CCDNUM":6 , "Offset": 0.0}
        infoDict["S28"] = {"xCenter": -152.172,"yCenter":  95.835, "FAflag":False, "CCDNUM":7 , "Offset": 0.0}
        infoDict["S29"] = {"xCenter": -185.988,"yCenter": -63.890, "FAflag":False, "CCDNUM":1 , "Offset": 0.0}
        infoDict["S30"] = {"xCenter": -185.988,"yCenter":   0.000, "FAflag":False, "CCDNUM":2 , "Offset": 0.0}
        infoDict["S31"] = {"xCenter": -185.988,"yCenter":  63.890, "FAflag":False, "CCDNUM":3 , "Offset": 0.0}
        infoDict["N1"]  = {"xCenter": 16.908,  "yCenter":-191.670, "FAflag":False, "CCDNUM":32, "Offset": 0.0}
        infoDict["N2"]  = {"xCenter": 16.908,  "yCenter":-127.780, "FAflag":False, "CCDNUM":33, "Offset": 0.0}
        infoDict["N3"]  = {"xCenter": 16.908,  "yCenter": -63.890, "FAflag":False, "CCDNUM":34, "Offset": 0.0}
        infoDict["N4"]  = {"xCenter": 16.908,  "yCenter":   0.000, "FAflag":False, "CCDNUM":35, "Offset": 0.0}
        infoDict["N5"]  = {"xCenter": 16.908,  "yCenter":  63.890, "FAflag":False, "CCDNUM":36, "Offset": 0.0}
        infoDict["N6"]  = {"xCenter": 16.908,  "yCenter": 127.780, "FAflag":False, "CCDNUM":37, "Offset": 0.0}
        infoDict["N7"]  = {"xCenter": 16.908,  "yCenter": 191.670, "FAflag":False, "CCDNUM":38, "Offset": 0.0}
        infoDict["N8"]  = {"xCenter": 50.724,  "yCenter":-159.725, "FAflag":False, "CCDNUM":39, "Offset": 0.0}
        infoDict["N9"]  = {"xCenter": 50.724,  "yCenter": -95.835, "FAflag":False, "CCDNUM":40, "Offset": 0.0}
        infoDict["N10"] = {"xCenter": 50.724,  "yCenter": -31.945, "FAflag":False, "CCDNUM":41, "Offset": 0.0}
        infoDict["N11"] = {"xCenter": 50.724,  "yCenter":  31.945, "FAflag":False, "CCDNUM":42, "Offset": 0.0}
        infoDict["N12"] = {"xCenter": 50.724,  "yCenter":  95.835, "FAflag":False, "CCDNUM":43, "Offset": 0.0}
        infoDict["N13"] = {"xCenter": 50.724,  "yCenter": 159.725, "FAflag":False, "CCDNUM":44, "Offset": 0.0}
        infoDict["N14"] = {"xCenter": 84.540,  "yCenter":-159.725, "FAflag":False, "CCDNUM":45, "Offset": 0.0}
        infoDict["N15"] = {"xCenter": 84.540,  "yCenter": -95.835, "FAflag":False, "CCDNUM":46, "Offset": 0.0}
        infoDict["N16"] = {"xCenter": 84.540,  "yCenter": -31.945, "FAflag":False, "CCDNUM":47, "Offset": 0.0}
        infoDict["N17"] = {"xCenter": 84.540,  "yCenter":  31.945, "FAflag":False, "CCDNUM":48, "Offset": 0.0}
        infoDict["N18"] = {"xCenter": 84.540,  "yCenter":  95.835, "FAflag":False, "CCDNUM":49, "Offset": 0.0}
        infoDict["N19"] = {"xCenter": 84.540,  "yCenter": 159.725, "FAflag":False, "CCDNUM":50, "Offset": 0.0}
        infoDict["N20"] = {"xCenter": 118.356, "yCenter":-127.780, "FAflag":False, "CCDNUM":51, "Offset": 0.0}
        infoDict["N21"] = {"xCenter": 118.356, "yCenter": -63.890, "FAflag":False, "CCDNUM":52, "Offset": 0.0}
        infoDict["N22"] = {"xCenter": 118.356, "yCenter":   0.000, "FAflag":False, "CCDNUM":53, "Offset": 0.0}
        infoDict["N23"] = {"xCenter": 118.356, "yCenter":  63.890, "FAflag":False, "CCDNUM":54, "Offset": 0.0}
        infoDict["N24"] = {"xCenter": 118.356, "yCenter": 127.780, "FAflag":False, "CCDNUM":55, "Offset": 0.0}
        infoDict["N25"] = {"xCenter": 152.172, "yCenter": -95.835, "FAflag":False, "CCDNUM":56, "Offset": 0.0}
        infoDict["N26"] = {"xCenter": 152.172, "yCenter": -31.945, "FAflag":False, "CCDNUM":57, "Offset": 0.0}
        infoDict["N27"] = {"xCenter": 152.172, "yCenter":  31.945, "FAflag":False, "CCDNUM":58, "Offset": 0.0}
        infoDict["N28"] = {"xCenter": 152.172, "yCenter":  95.835, "FAflag":False, "CCDNUM":59, "Offset": 0.0}
        infoDict["N29"] = {"xCenter": 185.988, "yCenter": -63.890, "FAflag":False, "CCDNUM":60, "Offset": 0.0}
        infoDict["N30"] = {"xCenter": 185.988, "yCenter":   0.000, "FAflag":False, "CCDNUM":61, "Offset": 0.0}
        infoDict["N31"] = {"xCenter": 185.988, "yCenter":  63.890, "FAflag":False, "CCDNUM":62, "Offset": 0.0}
        infoDict["FS1"] = {"xCenter": -152.172,"yCenter": 143.7525,"FAflag":True , "CCDNUM":66, "Offset": 1500.0}
        infoDict["FS2"] = {"xCenter": -185.988,"yCenter": 111.8075,"FAflag":True , "CCDNUM":65, "Offset":-1500.0}
        infoDict["FS3"] = {"xCenter": -219.804,"yCenter":  15.9725,"FAflag":True , "CCDNUM":63, "Offset": 1500.0}
        infoDict["FS4"] = {"xCenter": -219.804,"yCenter": -15.9725,"FAflag":True , "CCDNUM":64, "Offset":-1500.0}
        infoDict["FN1"] = {"xCenter": 152.172, "yCenter": 143.7525,"FAflag":True , "CCDNUM":67, "Offset":-1500.0}
        infoDict["FN2"] = {"xCenter": 185.988, "yCenter": 111.8075,"FAflag":True , "CCDNUM":68, "Offset": 1500.0}
        infoDict["FN3"] = {"xCenter": 219.804, "yCenter":  15.9725,"FAflag":True , "CCDNUM":69, "Offset":-1500.0}
        infoDict["FN4"] = {"xCenter": 219.804, "yCenter": -15.9725,"FAflag":True , "CCDNUM":70, "Offset": 1500.0}

        return infoDict
    

    def __init__(self,**inputDict):

        self.infoDict = self.info()
        self.mmperpixel = 0.015
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

        # CCD size in pixels
        if ccdinfo["FAflag"]:
            xpixHalfSize = 1024.
            ypixHalfSize = 1024.
        else:
            xpixHalfSize = 1024.
            ypixHalfSize = 2048.

        # calculate positions
        xPos = ccdinfo["xCenter"] + (float(ix)-xpixHalfSize+0.5)*self.mmperpixel
        yPos = ccdinfo["yCenter"] + (float(iy)-ypixHalfSize+0.5)*self.mmperpixel

        return xPos,yPos

    def getPixel(self,extname,xPos,yPos):
        """ given a coordinate in [mm], return pixel number
        """

        ccdinfo = self.infoDict[extname]

        # CCD size in pixels
        if ccdinfo["FAflag"]:
            xpixHalfSize = 1024.
            ypixHalfSize = 1024.
        else:
            xpixHalfSize = 1024.
            ypixHalfSize = 2048.

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
            if ccdinfo["FAflag"]:
                xpixHalfSize = 1024.
                ypixHalfSize = 1024.
            else:
                xpixHalfSize = 1024.
                ypixHalfSize = 2048.

            if nxdif <= xpixHalfSize and nydif <= ypixHalfSize:
                return ext

        # get to here if we are not inside a chip
        return None
            
        


class mosaicinfo(object):
    """ mosaicinfo is a class used to contain Mosaic geometry information and various utility routines
    """

    def info(self):
        # info returns a dictionary chock full of info on the DECam geometry
        # keyed by the CCD name

        infoDict = {}

        # store a dictionary for each CCD, keyed by the CCD name
        infoDict["0"] =  {"xCenter":  -(2048+1024)*self.mmperpixel,"yCenter":-2048*self.mmperpixel,"FAflag":False}
        infoDict["1"] =  {"xCenter":  -(1024)*self.mmperpixel,     "yCenter":-2048*self.mmperpixel,"FAflag":False}
        infoDict["2"] =  {"xCenter":   (1024)*self.mmperpixel,     "yCenter":-2048*self.mmperpixel,"FAflag":False}
        infoDict["3"] =  {"xCenter":   (2048+1024)*self.mmperpixel,"yCenter":-2048*self.mmperpixel,"FAflag":False}
        infoDict["4"] =  {"xCenter":  -(2048+1024)*self.mmperpixel,"yCenter": 2048*self.mmperpixel,"FAflag":False}
        infoDict["5"] =  {"xCenter":  -(1024)*self.mmperpixel,     "yCenter": 2048*self.mmperpixel,"FAflag":False}
        infoDict["6"] =  {"xCenter":   (1024)*self.mmperpixel,     "yCenter": 2048*self.mmperpixel,"FAflag":False}
        infoDict["7"] =  {"xCenter":   (2048+1024)*self.mmperpixel,"yCenter": 2048*self.mmperpixel,"FAflag":False}


        return infoDict
    

    def __init__(self,**inputDict):

        self.mmperpixel = 0.015
        self.infoDict = self.info()


    def getPosition(self,extname,ix,iy):
        # return the x,y position in [mm] for a given CCD and pixel number
        # note that the ix,iy are Image pixels - overscans removed - and start at zero

        ccdinfo = self.infoDict[extname]

        # CCD size in pixels
        xpixHalfSize = 1024.
        ypixHalfSize = 2048.

        # calculate positions
        xPos = ccdinfo["xCenter"] + (float(ix)-xpixHalfSize+0.5)*self.mmperpixel
        yPos = ccdinfo["yCenter"] + (float(iy)-ypixHalfSize+0.5)*self.mmperpixel

        return xPos,yPos
