#
# Utility methods for DECam
#
import numpy as np
from collections import OrderedDict

class decaminfo(object):
    """ decaminfo is a class used to contain DECam geometry information and various utility routines

    Aaron Roodman (C) SLAC National Accelerator Laboratory, Stanford University 2012.
    """

    def info(self):
        """info returns a dictionary chock full of info on the DECam geometry
        keyed by the CCD name
        """

        infoDict = OrderedDict()

        # store a dictionary for each CCD, keyed by the CCD name
        # AJR 9/14/2012 fixed these to agree with the DS9 coordinate system
        infoDict["S1"] =  {"xCenter":  -16.908,"yCenter":-191.670, "FAflag":False, "CCDNUM":25, "Offset": 0.0, "Extension": 1 }
        infoDict["S2"]  = {"xCenter":  -16.908,"yCenter":-127.780, "FAflag":False, "CCDNUM":26, "Offset": 0.0, "Extension": 2 }
        infoDict["S3"]  = {"xCenter":  -16.908,"yCenter": -63.890, "FAflag":False, "CCDNUM":27, "Offset": 0.0, "Extension": 3 } 
        infoDict["S4"]  = {"xCenter":  -16.908,"yCenter":   0.000, "FAflag":False, "CCDNUM":28, "Offset": 0.0, "Extension": 41 }
        infoDict["S5"]  = {"xCenter":  -16.908,"yCenter":  63.890, "FAflag":False, "CCDNUM":29, "Offset": 0.0, "Extension": 42 }
        infoDict["S6"]  = {"xCenter":  -16.908,"yCenter": 127.780, "FAflag":False, "CCDNUM":30, "Offset": 0.0, "Extension": 43 }
        infoDict["S7"]  = {"xCenter":  -16.908,"yCenter": 191.670, "FAflag":False, "CCDNUM":31, "Offset": 0.0, "Extension": 44 }
        infoDict["S8"]  = {"xCenter":  -50.724,"yCenter":-159.725, "FAflag":False, "CCDNUM":19, "Offset": 0.0, "Extension": 7 }
        infoDict["S9"]  = {"xCenter":  -50.724,"yCenter": -95.835, "FAflag":False, "CCDNUM":20, "Offset": 0.0, "Extension": 8 }
        infoDict["S10"] = {"xCenter":  -50.724,"yCenter": -31.945, "FAflag":False, "CCDNUM":21, "Offset": 0.0, "Extension": 19 }
        infoDict["S11"] = {"xCenter":  -50.724,"yCenter":  31.945, "FAflag":False, "CCDNUM":22, "Offset": 0.0, "Extension": 20 }
        infoDict["S12"] = {"xCenter":  -50.724,"yCenter":  95.835, "FAflag":False, "CCDNUM":23, "Offset": 0.0, "Extension": 21 }
        infoDict["S13"] = {"xCenter":  -50.724,"yCenter": 159.725, "FAflag":False, "CCDNUM":24, "Offset": 0.0, "Extension": 22 }
        infoDict["S14"] = {"xCenter":  -84.540,"yCenter":-159.725, "FAflag":False, "CCDNUM":13, "Offset": 0.0, "Extension": 9 }
        infoDict["S15"] = {"xCenter":  -84.540,"yCenter": -95.835, "FAflag":False, "CCDNUM":14, "Offset": 0.0, "Extension": 10 }
        infoDict["S16"] = {"xCenter":  -84.540,"yCenter": -31.945, "FAflag":False, "CCDNUM":15, "Offset": 0.0, "Extension": 25 }
        infoDict["S17"] = {"xCenter":  -84.540,"yCenter":  31.945, "FAflag":False, "CCDNUM":16, "Offset": 0.0, "Extension": 26 }
        infoDict["S18"] = {"xCenter":  -84.540,"yCenter":  95.835, "FAflag":False, "CCDNUM":17, "Offset": 0.0, "Extension": 23 }
        infoDict["S19"] = {"xCenter":  -84.540,"yCenter": 159.725, "FAflag":False, "CCDNUM":18, "Offset": 0.0, "Extension": 24 }
        infoDict["S20"] = {"xCenter": -118.356,"yCenter":-127.780, "FAflag":False, "CCDNUM":8 , "Offset": 0.0, "Extension": 11 }
        infoDict["S21"] = {"xCenter": -118.356,"yCenter": -63.890, "FAflag":False, "CCDNUM":9 , "Offset": 0.0, "Extension": 27 }
        infoDict["S22"] = {"xCenter": -118.356,"yCenter":   0.000, "FAflag":False, "CCDNUM":10, "Offset": 0.0, "Extension": 28 }
        infoDict["S23"] = {"xCenter": -118.356,"yCenter":  63.890, "FAflag":False, "CCDNUM":11, "Offset": 0.0, "Extension": 29 }
        infoDict["S24"] = {"xCenter": -118.356,"yCenter": 127.780, "FAflag":False, "CCDNUM":12, "Offset": 0.0, "Extension": 30 }
        infoDict["S25"] = {"xCenter": -152.172,"yCenter": -95.835, "FAflag":False, "CCDNUM":4 , "Offset": 0.0, "Extension": 12 }
        infoDict["S26"] = {"xCenter": -152.172,"yCenter": -31.945, "FAflag":False, "CCDNUM":5 , "Offset": 0.0, "Extension": 31 }
        infoDict["S27"] = {"xCenter": -152.172,"yCenter":  31.945, "FAflag":False, "CCDNUM":6 , "Offset": 0.0, "Extension": 32 }
        infoDict["S28"] = {"xCenter": -152.172,"yCenter":  95.835, "FAflag":False, "CCDNUM":7 , "Offset": 0.0, "Extension": 33 }
        infoDict["S29"] = {"xCenter": -185.988,"yCenter": -63.890, "FAflag":False, "CCDNUM":1 , "Offset": 0.0, "Extension": 34 }
        infoDict["S30"] = {"xCenter": -185.988,"yCenter":   0.000, "FAflag":False, "CCDNUM":2 , "Offset": 0.0, "Extension": 35 }
        infoDict["S31"] = {"xCenter": -185.988,"yCenter":  63.890, "FAflag":False, "CCDNUM":3 , "Offset": 0.0, "Extension": 36 }
        infoDict["N1"]  = {"xCenter": 16.908,  "yCenter":-191.670, "FAflag":False, "CCDNUM":32, "Offset": 0.0, "Extension": 4 }
        infoDict["N2"]  = {"xCenter": 16.908,  "yCenter":-127.780, "FAflag":False, "CCDNUM":33, "Offset": 0.0, "Extension": 5 }
        infoDict["N3"]  = {"xCenter": 16.908,  "yCenter": -63.890, "FAflag":False, "CCDNUM":34, "Offset": 0.0, "Extension": 6 }
        infoDict["N4"]  = {"xCenter": 16.908,  "yCenter":   0.000, "FAflag":False, "CCDNUM":35, "Offset": 0.0, "Extension": 37 }
        infoDict["N5"]  = {"xCenter": 16.908,  "yCenter":  63.890, "FAflag":False, "CCDNUM":36, "Offset": 0.0, "Extension": 38 }
        infoDict["N6"]  = {"xCenter": 16.908,  "yCenter": 127.780, "FAflag":False, "CCDNUM":37, "Offset": 0.0, "Extension": 39 }
        infoDict["N7"]  = {"xCenter": 16.908,  "yCenter": 191.670, "FAflag":False, "CCDNUM":38, "Offset": 0.0, "Extension": 40 }
        infoDict["N8"]  = {"xCenter": 50.724,  "yCenter":-159.725, "FAflag":False, "CCDNUM":39, "Offset": 0.0, "Extension": 13 }
        infoDict["N9"]  = {"xCenter": 50.724,  "yCenter": -95.835, "FAflag":False, "CCDNUM":40, "Offset": 0.0, "Extension": 14 }
        infoDict["N10"] = {"xCenter": 50.724,  "yCenter": -31.945, "FAflag":False, "CCDNUM":41, "Offset": 0.0, "Extension": 45 }
        infoDict["N11"] = {"xCenter": 50.724,  "yCenter":  31.945, "FAflag":False, "CCDNUM":42, "Offset": 0.0, "Extension": 46 }
        infoDict["N12"] = {"xCenter": 50.724,  "yCenter":  95.835, "FAflag":False, "CCDNUM":43, "Offset": 0.0, "Extension": 47 }
        infoDict["N13"] = {"xCenter": 50.724,  "yCenter": 159.725, "FAflag":False, "CCDNUM":44, "Offset": 0.0, "Extension": 48 }
        infoDict["N14"] = {"xCenter": 84.540,  "yCenter":-159.725, "FAflag":False, "CCDNUM":45, "Offset": 0.0, "Extension": 15 }
        infoDict["N15"] = {"xCenter": 84.540,  "yCenter": -95.835, "FAflag":False, "CCDNUM":46, "Offset": 0.0, "Extension": 16 }
        infoDict["N16"] = {"xCenter": 84.540,  "yCenter": -31.945, "FAflag":False, "CCDNUM":47, "Offset": 0.0, "Extension": 51 }
        infoDict["N17"] = {"xCenter": 84.540,  "yCenter":  31.945, "FAflag":False, "CCDNUM":48, "Offset": 0.0, "Extension": 52 }
        infoDict["N18"] = {"xCenter": 84.540,  "yCenter":  95.835, "FAflag":False, "CCDNUM":49, "Offset": 0.0, "Extension": 49 }
        infoDict["N19"] = {"xCenter": 84.540,  "yCenter": 159.725, "FAflag":False, "CCDNUM":50, "Offset": 0.0, "Extension": 50 }
        infoDict["N20"] = {"xCenter": 118.356, "yCenter":-127.780, "FAflag":False, "CCDNUM":51, "Offset": 0.0, "Extension": 17 }
        infoDict["N21"] = {"xCenter": 118.356, "yCenter": -63.890, "FAflag":False, "CCDNUM":52, "Offset": 0.0, "Extension": 53 }
        infoDict["N22"] = {"xCenter": 118.356, "yCenter":   0.000, "FAflag":False, "CCDNUM":53, "Offset": 0.0, "Extension": 54 }
        infoDict["N23"] = {"xCenter": 118.356, "yCenter":  63.890, "FAflag":False, "CCDNUM":54, "Offset": 0.0, "Extension": 55 }
        infoDict["N24"] = {"xCenter": 118.356, "yCenter": 127.780, "FAflag":False, "CCDNUM":55, "Offset": 0.0, "Extension": 56 }
        infoDict["N25"] = {"xCenter": 152.172, "yCenter": -95.835, "FAflag":False, "CCDNUM":56, "Offset": 0.0, "Extension": 18 }
        infoDict["N26"] = {"xCenter": 152.172, "yCenter": -31.945, "FAflag":False, "CCDNUM":57, "Offset": 0.0, "Extension": 57 }
        infoDict["N27"] = {"xCenter": 152.172, "yCenter":  31.945, "FAflag":False, "CCDNUM":58, "Offset": 0.0, "Extension": 58 }
        infoDict["N28"] = {"xCenter": 152.172, "yCenter":  95.835, "FAflag":False, "CCDNUM":59, "Offset": 0.0, "Extension": 59 }
        infoDict["N29"] = {"xCenter": 185.988, "yCenter": -63.890, "FAflag":False, "CCDNUM":60, "Offset": 0.0, "Extension": 60 }
        infoDict["N30"] = {"xCenter": 185.988, "yCenter":   0.000, "FAflag":False, "CCDNUM":61, "Offset": 0.0, "Extension": 61 }
        infoDict["N31"] = {"xCenter": 185.988, "yCenter":  63.890, "FAflag":False, "CCDNUM":62, "Offset": 0.0, "Extension": 62 }
        infoDict["FS1"] = {"xCenter": -152.172,"yCenter": 143.7525,"FAflag":True , "CCDNUM":66, "Offset": 1500.0, "Extension": 63 }
        infoDict["FS2"] = {"xCenter": -185.988,"yCenter": 111.8075,"FAflag":True , "CCDNUM":65, "Offset":-1500.0, "Extension": 64 }
        infoDict["FS3"] = {"xCenter": -219.804,"yCenter":  15.9725,"FAflag":True , "CCDNUM":63, "Offset": 1500.0, "Extension": 65 }
        infoDict["FS4"] = {"xCenter": -219.804,"yCenter": -15.9725,"FAflag":True , "CCDNUM":64, "Offset":-1500.0, "Extension": 66 }
        infoDict["FN1"] = {"xCenter": 152.172, "yCenter": 143.7525,"FAflag":True , "CCDNUM":67, "Offset":-1500.0, "Extension": 67 }
        infoDict["FN2"] = {"xCenter": 185.988, "yCenter": 111.8075,"FAflag":True , "CCDNUM":68, "Offset": 1500.0, "Extension": 68 }
        infoDict["FN3"] = {"xCenter": 219.804, "yCenter":  15.9725,"FAflag":True , "CCDNUM":69, "Offset":-1500.0, "Extension": 69 }
        infoDict["FN4"] = {"xCenter": 219.804, "yCenter": -15.9725,"FAflag":True , "CCDNUM":70, "Offset": 1500.0, "Extension": 70 }

        return infoDict
    
    def mkArraybyNum(self):
        infoArr = np.zeros((71,3)) # so that x,y,faflag = infoArr[ccdnum,:] 
        for ccdname in self.infoDict.keys():
            infoArr[self.infoDict[ccdname]['CCDNUM'],0] = self.infoDict[ccdname]['xCenter']
            infoArr[self.infoDict[ccdname]['CCDNUM'],1] = self.infoDict[ccdname]['yCenter']
            if self.infoDict[ccdname]['FAflag']:
                infoArr[self.infoDict[ccdname]['CCDNUM'],2] = 1
            else:
                infoArr[self.infoDict[ccdname]['CCDNUM'],2] = 0

        return infoArr

    def __init__(self,**inputDict):

        self.infoDict = self.info()
        self.mmperpixel = 0.015
        self.rClear = 225.0
        self.infoArrbyNum = self.mkArraybyNum()

    def __getstate__(self):
        stateDict = {}
        keysToPickle = ['infoDict','mmperpixel','rClear','infoArrbyNum']
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

    def getPositionByNum(self,ccdnum,ix,iy):
        """ return the x,y position in [mm] for a given CCD and pixel number                                                                                       
        note that the ix,iy are Image pixels - overscans removed - and start at zero                                                                            
        make sure this code is vectorizable...
        """

        # CCD size in pixels, special code for Focus sensors
        xpixHalfSize = 1024.
        ypixHalfSize = np.where(self.infoArrbyNum[ccdnum,2]==1,1024.,2048.)

        # calculate positions                                                                                                                                      
        xPos = self.infoArrbyNum[ccdnum,0] + (float(ix)-xpixHalfSize+0.5)*self.mmperpixel
        yPos = self.infoArrbyNum[ccdnum,1] + (float(iy)-ypixHalfSize+0.5)*self.mmperpixel

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
        for ext in list(self.infoDict.keys()):
            ccdinfo = self.infoDict[ext]
            # is this x,y inside this chip?
            nxdif = np.abs( (xPos - ccdinfo["xCenter"]) / self.mmperpixel )
            nydif = np.abs( (yPos - ccdinfo["yCenter"]) / self.mmperpixel )

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
