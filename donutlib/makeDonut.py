#!/usr/bin/env python
# $Rev:: 203                                                          $:  
# $Author:: roodman                                                   $:  
# $LastChangedDate:: 2015-05-20 10:20:01 -0700 (Wed, 20 May 2015)     $: 
#
# Make one calculated image, either from Zemax Zernike array
# or from Zemax WFM file.  
#
import numpy
import scipy
from astropy.io import fits as pyfits
from array import array
import argparse
import pdb
import sys
import string
from donutlib.donututil import getZemaxWfm,getFitsWfm
from sandbox.donutengine import donutengine
from donutlib.decamutil import decaminfo

class makeDonut(object):
    """
    makeDonut is a class used make artificial donuts or stars
        
    Aaron Roodman - SLAC National Accelerator Laboratory
    """


    def updateParam(**inputDict):
        """ set parameters
        """
        paramDict = {"inputFile":"",
                 "wfmFile":"",
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
                 "gain":1.0,
                 "flipFlag":False,
                 "ZernikeArray":[]}

        paramDict.update(inputDict)
        return paramDict
    
    def __init__(self,**inputDict):
        """ initialize and build donut maker
        """
        paramDict = updateParam(**inputDict)

        # check parameters are ok
        if paramDict["nbin"] != paramDict["nPixels"]*paramDict["pixelOverSample"]:
            print "makeDonut:  nbin must = nPixels * pixelOverSample !!!"
            sys.exit(1)

        # also require that _Lu > 2R, ie. that pupil fits in pupil plane
        # this translates into requiring that (Lambda F / pixelSize) * pixelOverSample * scaleFactor > 1
        # Why does this depend on scaleFactor, but not Z4?? Answer: scaleFactor effectively changes the wavelength, so 
        # it must be included.  It is also possible that Z4 is too big for the nPixels - buts that another limit than this one
        F = 2.9  # hardcode for DECam for now
        pixelSize = 15.e-6 
        if paramDict["pixelOverSample"] * paramDict["scaleFactor"] * (paramDict["waveLength"] * F / pixelSize) < 1. :
            print "makeDonut:  ERROR pupil doesn't fit!!!"
            print "            value = ",paramDict["pixelOverSample"] * paramDict["scaleFactor"] * (paramDict["waveLength"] * F / pixelSize)
            sys.exit(2)

        # for WFM, need to turn gridCalcMode to False for donutengine
        #if paramDict["useWFM"]:
        #    paramDict["gridCalcMode"] = False
        #
        # don't do this by default, may need it if Zemax is used, but not if correct bins sizes are used.
        # and it isn't coded with c++ donutengine yet either...

        # declare fit function
        self.gFitFunc = donutengine(**paramDict)


    def generate(**inputDict):
        """ generate a donut or star
        """
        paramDict = updateParam(**inputDict)
        
        # call setXYDECam(x,y) to get an off-axis pupil function!!!
        self.gFitFunc.setXYDECam(paramDict["xDECam"],paramDict["yDECam"])

        # convert this position to extname,ix,iy
        # Note: x=0,y=0 is in between sensors.  IF this is used, then just set these by hand
        if paramDict["iTelescope"]==0:
            dinfo = decaminfo()
            if paramDict["xDECam"]==0.0 and paramDict["yDECam"]==0.0:
                extname = "N4"
                ix = -1024
                iy = 2048
            else:        
                extname = dinfo.getSensor(paramDict["xDECam"],paramDict["yDECam"])
                x,y = dinfo.getPixel(extname,paramDict["xDECam"],paramDict["yDECam"])
                ix = int(x+0.5)
                iy = int(y+0.5)
        else:
            extname = "Other"
            ix = 0
            iy = 0

        # parameters for Donuts
        par = numpy.zeros(self.gFitFunc.npar)
        par[self.gFitFunc.ipar_rzero] = paramDict["rzero"]
        par[self.gFitFunc.ipar_nEle] = paramDict["nEle"]
        par[self.gFitFunc.ipar_bkgd] = paramDict["background"]

        # get command line input for Zernikes (if any)
        inputZernikeArray  = numpy.array(paramDict["ZernikeArray"])
        someInputZernike = inputZernikeArray.any()

        # WFM or Zernikes
        if paramDict["wfmFile"]=="" or someInputZernike :

            # test for a Zernike input file or input from cmd line, or both
            if paramDict["inputFile"] != ""  :

                aDir,aFile = os.path.split(paramDict["inputFile"])
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

        elif paramDict["wfmFile"]!=""  :

            # get WFM array from the file
            if string.find(paramDict["wfmFile"],".txt") >= 0:
                xaxis,yaxis,wfm = getZemaxWfm(paramDict["wfmFile"])
            elif string.find(paramDict["wfmFile"],".fits") >= 0:
                wfm = getFitsWfm(paramDict["wfmFile"],0)
                # now make the donut
                self.gFitFunc.fillPar(par)
                self.gFitFunc.calcWFMtoImage(wfm)

        # did we get this far?
        #print "makeDonut: calcAll is finished!"


        # randomize
        theImage = self.gFitFunc.getvImage().copy()
        if paramDict["randomFlag"]:
            numpy.random.seed(paramDict["randomSeed"])
            postageshape = theImage.shape
            nranval = numpy.random.normal(0.0,1.0,postageshape)
            imarr = theImage + nranval*numpy.sqrt(theImage)
        else:
            imarr = theImage

        # apply gain
        imarr = imarr / paramDict["gain"]

        # make sure that imarr has dtype = float32
        imarr = numpy.float32(imarr)

        # if desired flip array in X
        if paramDict["flipFlag"]:
            imarr = numpy.fliplr(imarr)

        # get Zernike code
        #myZern = self.gFitFunc.ZernikeObject.ZernikeDescription

        #output the calculated image

        #print "makeDonut: ready to make output file"

        # calculated Donut
        if paramDict["writeToFits"]:
            hdu = pyfits.PrimaryHDU(imarr)
            prihdr =  hdu.header

            prihdr.set("SCALE",F,"Arsec/pixel")
            prihdr.set("XDECAM",paramDict["xDECam"],"Target xposition (mm) in focal plane")
            prihdr.set("YDECAM",paramDict["yDECam"],"Target yposition (mm) in focal plane")
            prihdr.set("EXTNAME",extname)  
            prihdr.set("IX",ix)
            prihdr.set("IY",iy)
            prihdr.set("FILTER",3,"Filter number 1-6=ugrizY")
            prihdr.set("FILTNAME","r","Filter name")

            prihdr.set("GAIN",paramDict["gain"],"Np.e./ADU")
            prihdr.set("nEleInp",par[self.gFitFunc.ipar_nEle],"Number of photo-electrons")
            prihdr.set("rzeroInp",par[self.gFitFunc.ipar_rzero],"Fried Parameters [m]")
            prihdr.set("bkgdInp",par[self.gFitFunc.ipar_bkgd],"Background")

            if paramDict["wfmFile"]=="" :
                for iZ in range(self.gFitFunc.nZernikeSize):
                    name = "Z" + str(iZ+2)
                    prihdr.set(name,par[self.gFitFunc.ipar_ZernikeFirst+iZ])

            hdulist = pyfits.HDUList([hdu])
            outFile = paramDict["outputPrefix"] + ".stamp.fits"
            hdulist.writeto(outFile,clobber=True)

        else:
            return imarr



        #print "makeDonut is DONE!"

#
#  if running from the command line, then we need the ArgumentParser, otherwise
#  bring in the module makeDonut
#
if __name__ == "__main__":

    parser = argparse.ArgumentParser(prog='makeDonut')
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
    m = makeDonut(**aDict)
    m.generate(**aDict)