#!/usr/bin/env python
# $Rev:: 212                                                          $:  
# $Author:: roodman                                                   $:  
# $LastChangedDate:: 2015-08-17 12:07:18 -0700 (Mon, 17 Aug 2015)     $: 
#
# Make one calculated image, either from Zemax Zernike array
# or from Zemax WFM file.  
#
import numpy
import scipy
#import pyfits
from astropy.io import fits as pyfits
from array import array
import argparse
import pdb
import sys
import string
from donutlib.donututil import getZemaxWfm,getFitsWfm
#from donutengine import donutengine
from donutlib.donutengine import donutengine
from donutlib.decamutil import decaminfo


class makedonut(object):
    """ make simulated Donuts

    Aaron Roodman (C) SLAC National Accelerator Laboratory, Stanford University 2012.
    """

    def __init__(self,**inputDict):
        """ initialize
        """

        # initialize the parameter Dictionary, and update defaults from inputDict
        self.paramDict = {"inputFile":"",
                          "wfmFile":"",
                          "wfmArray":None,
                          "writeToFits":False,
                          "outputPrefix":"testone",
                          "xDECam":0.0,
                          "yDECam":0.0,
                          "debugFlag":False,
                          "rootFlag":False,
                          "iTelescope":0,
                          "waveLength":700.0e-9,
                          "nZernikeTerms":37,
                          "nbin":512,
                          "nPixels":64,
                          "gridCalcMode":True,
                          "pixelOverSample":8,
                          "scaleFactor":1.,                 
                          "rzero":0.125,
                          "nEle":1.0e6,
                          "background":4000.,
                          "randomFlag":False,
                          "randomSeed":209823,  # if this is an invalid integer, crazy errors will ensue
                          "gain":1.0,
                          "flipFlag":False,
                          "ZernikeArray":[]}

        self.paramDict.update(inputDict)

        # check parameters are ok
        if self.paramDict["nbin"] != self.paramDict["nPixels"]*self.paramDict["pixelOverSample"]:
            print "makedonut:  nbin must = nPixels * pixelOverSample !!!"
            sys.exit(1)

        # also require that _Lu > 2R, ie. that pupil fits in pupil plane
        # this translates into requiring that (Lambda F / pixelSize) * pixelOverSample * scaleFactor > 1
        # Why does this depend on scaleFactor, but not Z4?? Answer: scaleFactor effectively changes the wavelength, so 
        # it must be included.  It is also possible that Z4 is too big for the nPixels - buts that another limit than this one
        F = 2.9  # hardcode for DECam for now
        pixelSize = 15.e-6 
        if self.paramDict["pixelOverSample"] * self.paramDict["scaleFactor"] * (self.paramDict["waveLength"] * F / pixelSize) < 1. :
            print "makedonut:  ERROR pupil doesn't fit!!!"
            print "            value = ",self.paramDict["pixelOverSample"] * self.paramDict["scaleFactor"] * (self.paramDict["waveLength"] * F / pixelSize)
            #sys.exit(2)

        # for WFM, need to turn gridCalcMode to False for donutengine
        #if self.paramDict["useWFM"]:
        #    self.paramDict["gridCalcMode"] = False
            #
            # don't do this by default, may need it if Zemax is used, but not if correct bins sizes are used.
            # and it isn't coded with c++ donutengine yet either...
            
        # declare fit function
        self.gFitFunc = donutengine(**self.paramDict)

        # DECam info
        self.dinfo = decaminfo()

        # set random seed
        numpy.random.seed(self.paramDict["randomSeed"])


    def setXY(self,X,Y):
        # call setXYDECam(x,y) to get an off-axis pupil function!!!
        self.gFitFunc.setXYDECam(X,Y)

        # convert this position to extname,ix,iy
        # Note: x=0,y=0 is in between sensors.  IF this is used, then just set these by hand
        if X==0.0 and Y==0.0:
            self.extname = "N4"
            self.ix = -1024
            self.iy = 2048
        else:        
            self.extname = dinfo.getSensor(X,Y)
            x,y = dinfo.getPixel(self.extname,X,Y)
            self.ix = int(x+0.5)
            self.iy = int(y+0.5)

    def make(self,**inputDict):
        """ make the Donut
        """

        # update the parameters
        self.paramDict.update(inputDict)

        # update the focal plane position
        self.setXY(self.paramDict["xDECam"],self.paramDict["yDECam"])

        # parameters for Donuts
        par = numpy.zeros(self.gFitFunc.npar)
        par[self.gFitFunc.ipar_rzero] = self.paramDict["rzero"]
        par[self.gFitFunc.ipar_nEle] = self.paramDict["nEle"]
        par[self.gFitFunc.ipar_bkgd] = self.paramDict["background"]

        # get command line input for Zernikes (if any)
        inputZernikeArray  = numpy.array(self.paramDict["ZernikeArray"])
        someInputZernike = inputZernikeArray.any()

        # WFM or Zernikes
        if self.paramDict["wfmFile"]=="" or someInputZernike :

            # test for a Zernike input file or input from cmd line, or both
            if self.paramDict["inputFile"] != ""  :

                aDir,aFile = os.path.split(self.paramDict["inputFile"])
                if not sys.path.__contains__(aDir):
                    sys.path.append(aDir)

                # check .py in file name
                if aFile[-3:] == ".py" :
                    aaFile = aFile[0:-3]
                else :
                    aaFile = aFile

                ZernCoeff = __import__(aaFile).ZernCoeff
                print ZernCoeff

                for iZ in range(self.gFitFunc.nZernikeSize):
                    par[self.gFitFunc.ipar_ZernikeFirst+iZ] = ZernCoeff[iZ+1]

            # cmd line inputs overwrite what was in the file  -- start from Zernike2 in the input 
            if someInputZernike:
                for iZ in range(len(inputZernikeArray)-1):            
                    par[self.gFitFunc.ipar_ZernikeFirst+iZ] += inputZernikeArray[iZ+1] 

            # now make the Donut
            self.gFitFunc.calcAll(par)

        elif self.paramDict["wfmFile"]!=""  :

            # get WFM array from the file
            if string.find(self.paramDict["wfmFile"],".txt")>=0:
                xaxis,yaxis,wfm = getZemaxWfm(self.paramDict["wfmFile"])
            elif string.find(self.paramDict["wfmFile"],".fits")>=0:
                wfm = getFitsWfm(self.paramDict["wfmFile"],0)
                # now make the donut
                self.gFitFunc.fillPar(par)
                self.gFitFunc.calcWFMtoImage(wfm)

        elif self.paramDict["wfmArray"] != None:
            self.gFitFunc.fillPar(par)
            self.gFitFunc.calcWFMtoImage(self.paramDict["wfmArray"])

        # did we get this far?
        #print "makedonut: calcAll is finished!"

        # randomize
        theImage = self.gFitFunc.getvImage().copy()

        #print "makedonut: got the vImage"

        if self.paramDict["randomFlag"]:
            postageshape = theImage.shape
            nranval = numpy.random.normal(0.0,1.0,postageshape)
            imarr = theImage + nranval*numpy.sqrt(theImage)
        else:
            imarr = theImage

        # apply gain
        imarr = imarr / self.paramDict["gain"]

        # make sure that imarr has dtype = float32
        if self.paramDict["writeToFits"]:
            imarr = numpy.float32(imarr)
        else:
            imarr = numpy.float64(imarr)

        # if desired flip array in X
        if self.paramDict["flipFlag"]:
            imarr = numpy.fliplr(imarr)

        # get Zernike code
        #myZern = self.gFitFunc.ZernikeObject.ZernikeDescription

        #output the calculated image

        #print "makedonut: ready to make output file"

        # calculated Donut
        if self.paramDict["writeToFits"]:
            hdu = pyfits.PrimaryHDU(imarr)
            prihdr =  hdu.header

            prihdr.set("SCALE",0.27,"Arsec/pixel")
            prihdr.set("XDECAM",self.paramDict["xDECam"],"Target xposition (mm) in focal plane")
            prihdr.set("YDECAM",self.paramDict["yDECam"],"Target yposition (mm) in focal plane")
            prihdr.set("EXTNAME",self.extname)  
            prihdr.set("IX",self.ix)
            prihdr.set("IY",self.iy)
            prihdr.set("FILTER",3,"Filter number 1-6=ugrizY")
            prihdr.set("FILTNAME","r","Filter name")

            prihdr.set("GAIN",self.paramDict["gain"],"Np.e./ADU")
            prihdr.set("nEleInp",par[self.gFitFunc.ipar_nEle],"Number of photo-electrons")
            prihdr.set("rzeroInp",par[self.gFitFunc.ipar_rzero],"Fried Parameters [m]")
            prihdr.set("bkgdInp",par[self.gFitFunc.ipar_bkgd],"Background")

            if self.paramDict["wfmFile"]=="" :
                for iZ in range(self.gFitFunc.nZernikeSize):
                    name = "Z" + str(iZ+2)
                    prihdr.set(name,par[self.gFitFunc.ipar_ZernikeFirst+iZ])

            hdulist = pyfits.HDUList([hdu])
            outFile = self.paramDict["outputPrefix"] + ".stamp.fits"
            hdulist.writeto(outFile,clobber=True)
            return 1

        else:
            return imarr



#
#  if running from the command line, then we need the ArgumentParser, otherwise
#  bring in the module makeDonut
#
if __name__ == "__main__":

    parser = argparse.ArgumentParser(prog='makedonut')
    parser.add_argument("-i", "--inputFile",
                  dest="inputFile",
                  default="",
                  help="input file name")
    parser.add_argument("-w", "--wfmFile",
                  dest="wfmFile",
                  default="",
                  help="wfm input file name, .txt for Zemax, .fits or .fits.fz for Fits")
    parser.add_argument("-f", "--writeToFits",
                  dest="writeToFits",
                  action="store_true",
                  default=False,
                  help="set -f to write file out to a fits file, not the default")
    parser.add_argument("-o", "--outputPrefix",
                  dest="outputPrefix",
                  default="test",
                  help="output file prefix")
    parser.add_argument("-x", "--xDECam",
                  dest="xDECam",
                  default=0.,type=float,
                  help="focal plane x coordinate")
    parser.add_argument("-y", "--yDECam",
                  dest="yDECam",
                  default=0.,type=float,
                  help="focal plane y coordinate")    
    parser.add_argument("-d", "--debugFlag",
                  dest="debugFlag",
                  action="store_true",
                  default=False,
                  help="make debug plots too")
    parser.add_argument("-pl", "--printLevel",
                  dest="printLevel",
                  default=0,type=int,
                  help="debug printout")
    parser.add_argument("-root", "--rootFlag",
                  dest="rootFlag",
                  action="store_true",
                  default=False,
                  help="make root plots too")
    parser.add_argument("-flip", "--flipFlag",
                  dest="flipFlag",
                  action="store_true",
                  default=False,
                  help="if true flip the image in X, to mimic Zemax -> DECam")
    parser.add_argument("-t", "--iTelescope",
                  dest="iTelescope",
                  default=0,type=int,
                  help="iTelescope flag for DonutEngine, =0 DECam, =1 MosaicII")
    parser.add_argument("-wave", "--waveLength",
                  dest="waveLength",
                  default=700.0e-9,type=float,
                  help="waveLength [meters]")
    parser.add_argument("-nZ",
                dest="nZernikeTerms",
                default=37,type=int,
                help="number of Zernike terms")
    
    parser.add_argument("-nbin", "--nbin",
                  dest="nbin",
                  default=512,type=int,
                  help="Number of Bins")
    parser.add_argument("-nP", "--nPixels",
                  dest="nPixels",
                  default=64,type=int,
                  help="Number of Pixels")
    parser.add_argument("-gcmF", "--gridCalcMode",
                  dest="gridCalcMode",
                  action="store_false",
                  default=True,
                  help="turn off grid Calc mode in donutengine")
    parser.add_argument("-pos", "--pixelOverSample",
                  dest="pixelOverSample",
                  default=8,type=int,
                  help="pixel oversample factor")
    parser.add_argument("-s", "--scaleFactor",
                  dest="scaleFactor",
                  default=1.0,type=float,
                  help="scale Factor for FFT grid")

    parser.add_argument("-r","--rzero",
                    dest="rzero",
                    default=0.125,type=float,
                    help="Fried parameter")
    parser.add_argument("-n","--nEle",
                    dest="nEle",
                    default=1.0e6,type=float,
                    help="Number of photon-electrons")
    parser.add_argument("-b","--background",
                    dest="background",
                    default=4000.00,type=float,
                    help="background level")
    parser.add_argument("-ran","--randomFlag",
                    dest="randomFlag",
                    default=False,action="store_true",
                    help="randomize")
    parser.add_argument("-seed","--randomSeed",
                    dest="randomSeed",
                    default=2390487,type=int,
                    help="randomize")
    parser.add_argument("-g","--gain",
                    dest="gain",
                    default=1.00,type=float,
                    help="gain in Nphoto-electron/ADU")

    parser.add_argument("-za",
                    dest="ZernikeArray",
                    default=None,
                    help="""input Zernike Array with format "[Z1,Z2,Z3,Z4,Z5,Z6...Z37]" """)

    # collect the options
    options = parser.parse_args()
    aDict = vars(options)  #converts the object options to dictionary of key:value

    # need to convert ZernikeArray from text to an array, and also
    if aDict["ZernikeArray"]!=None:
        actualZernikeArray  = numpy.array(eval(aDict["ZernikeArray"]))
        aDict["ZernikeArray"] = actualZernikeArray
        print "Zernikes are", actualZernikeArray
    else:
        aDict["ZernikeArray"] = numpy.array([])
    print aDict

    # do it now
    m = makedonut(**aDict)
    m.make(**aDict)
