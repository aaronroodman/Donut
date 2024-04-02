import time
import numpy as np
from collections import OrderedDict
from astropy.io import fits as pyfits
import tabulate as tab
import pdb

from iminuit import Minuit
from donutlib.donutengine import donutengine
from donutlib.donututil import loadImage
from donutlib.donututil import calcStarting
from donutlib.decamutil import decaminfo

class donutfit(object):
    """ donutfit is a class used to fit donuts, using donutengine and MINUIT, for the DES experiment

    Aaron Roodman (C) SLAC National Accelerator Laboratory, Stanford University 2012.
    """

    def __init__(self,**inputDict):
        # init contains all initializations which are done only once for all fits
        # parameters in fixParamArray1 are nEle,rzero,bkgd,Z2,Z3,Z4,....Z11

        self.paramDict = {"nZernikeTerms":11,
                          "nFits":1,
                          "fixedParamArray1":[0,1,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1],  #need defaults up to quadrefoil
                          "debugFlag":False,
                          "outputWavefront":False,
                          "outputDiff":True,
                          "outputChi2":False,
                          "printLevel":1,
                          "maxIterations":1000,
                          "calcRzeroDerivative":True,
                          "wavefrontMap":None,
                          "deltaWFM":None,
                          "waveLength":700.0e-9,    # be sure this is in meters!
                          "doGridFit":False,
                          "spacing":64,
                          "gain":1.0}

        # search for key in inputDict, change defaults
        self.paramDict.update(inputDict)

        # setup the fit engine
        self.gFitFunc = donutengine(**self.paramDict)

        # need dummy versions before calling self.chisq
        self.imgarray = np.zeros(1)
        self.weight = np.ones(1)
        self.sigmasq = np.ones(1)

        # get decam info
        self.decamInfo = decaminfo()

        # do initial setup of Minuit parameters

        # status/limit arrays for Minuit parameters
        self.startingParam = np.zeros(self.gFitFunc.npar)
        self.errorParam = np.ones(self.gFitFunc.npar)
        self.loParam = np.zeros(self.gFitFunc.npar)
        self.hiParam = np.zeros(self.gFitFunc.npar)
        self.paramStatusArray = np.zeros(self.gFitFunc.npar)   # store =0 Floating, =1 Fixed

        # setup MINUIT
        self.gMinuit = Minuit(self.chisq,self.startingParam,name=self.gFitFunc.parNames,grad=self.calcgrad)
        self.gMinuit.strategy = 0
        self.gMinuit.errordef = Minuit.LEAST_SQUARES

        # set printlevel
        self.gMinuit.print_level = self.paramDict["printLevel"]

        # get wavefrontMap object
        self.wavefrontMap = self.paramDict['wavefrontMap']

        # if we have a deltaWFM, read it and install it
        if self.paramDict['deltaWFM'] and self.paramDict['deltaWFM']!='None':
            wfm_hdu = pyfits.open(self.paramDict['deltaWFM'])
            wfm = np.array(wfm_hdu[0].data,dtype=np.float64)  #need to convert when reading from fits!
            self.gFitFunc.setDeltaWFM(wfm)


    def chisq(self,par):

        # call donutengine to calculate image
        self.gFitFunc.calcAll(par)

        # compare to calculated image
        diff = self.imgarray - self.gFitFunc.getvImage()
        self.pullsq = diff*diff/self.sigmasq
        chisquared = self.pullsq.sum()

        # printout
        if self.paramDict["printLevel"]>=2:
            print('donutfit: Chi2 = ',chisquared)

        # save parameters for next iteration
        self.gFitFunc.savePar()

        # return result
        return chisquared

    def calcgrad(self,par):

        # make sure we've called chisq
        chisquared = self.chisq(par)

        # not currently called in calcAll
        self.gFitFunc.calcDerivatives(self.imgarray,self.weight)
        dChi2dpar = self.gFitFunc.getDerivatives()

        return dChi2dpar


    def setupFit(self,**inputFitDict):
        """ setup the fit, and do the fit, for a new Donut
        """

        # these defaults come from Zemax
        # and we use the entry keyed by None as a default if this is not from DECam
        # Fixed  FN's Z4 signs on 10/4/2012  AJR
        # still need to check the signs of the Trefoil terms...
        # Changed default zern4 to be 11.0, since fits prefer to start high 10/5/2012 AJR
        #
        inputZernikeDict = {'None':[0.0,0.0,0.0],
                            "FS1":[0.0,0.0,11.0,0.0,0.0,0.0,0.0,0.20,-0.17,-0.08],
                            "FS2":[0.0,0.0,-11.0,0.0,0.0,0.0,0.0,0.26,-0.01,-0.13],
                            "FS3":[0.0,0.0,11.0,0.0,0.0,0.0,0.0,0.05,0.25,-0.11],
                            "FS4":[0.0,0.0,-11.0,0.0,0.0,0.0,0.0,-0.05,0.25,-0.14],
                            "FN1":[0.0,0.0,-11.0,0.0,0.0,0.0,0.0,0.20,0.17,-0.08],
                            "FN2":[0.0,0.0,11.0,0.0,0.0,0.0,0.0,0.26,0.01,-0.13],
                            "FN3":[0.0,0.0,-11.0,0.0,0.0,0.0,0.0,0.05,-0.25,-0.11],
                            "FN4":[0.0,0.0,11.0,0.0,0.0,0.0,0.0,-0.05,-0.25,-0.14] }
        # default Dictionary
        self.fitDict = {"inputImageArray":None,
                        "inputFile":"",
                        "outputPrefix":"test",
                        "inputZernikeDict":inputZernikeDict,
                        "inputrzero":None,
                        "inputnEle":None,
                        "inputbkgd":None,
                        "inputZ4max":100.0,    #Sept 22, 2014 change default back to 100
                        "skipGradient":False}

        # update Dictionary with inputs
        self.fitDict.update(inputFitDict)

        # if desired print dict
        if self.paramDict["printLevel"]>=2:
            print(self.fitDict)

        # load image either from a file or from an input array
        if self.fitDict["inputImageArray"] ==None:

            # check for fits in file name
            fileName = self.fitDict["inputFile"]
            if fileName[-5:] != ".fits" :
                fileName = fileName + ".fits"

            # get the input file header
            hdulist = pyfits.open(fileName)
            iextension = 0   # set to other values for test usage
            self.inputHeader = hdulist[iextension].header

            # get the extension name from the header
            if list(self.inputHeader.keys()).count("EXTNAME")>0:
                extname = self.inputHeader["EXTNAME"]
                if extname == "":
                    extname = 'None'
            else:
                extname = 'None'

            # also get the IX,IY values from the header
            if list(self.inputHeader.keys()).count("IX")>0:
                ix = self.inputHeader["IX"]
            else:
                ix = 0.

            if list(self.inputHeader.keys()).count("IY")>0:
                iy = self.inputHeader["IY"]
            else:
                iy = 0.

            # calculate position in DECam focal plane
            # if this info isn't in the header, just set values to 0,0
            if extname != 'None':
                try:
                    xDECam,yDECam = self.decamInfo.getPosition(extname,ix,iy)
                except:
                    xDECam = 0.
                    yDECam = 0.
            else:
                xDECam = 0.
                yDECam = 0.

            # load the Image AJR 9/14/2012 - now assume we are only running on postage stamps - remove ability to
            # work on full image, never use that anymore...
            ### self.imgarray = hdulist[0].data.copy() this caused bugs with some fits files
            gain = self.paramDict["gain"] #convert to Nele from Nadu
            self.imgarray = gain * hdulist[iextension].data.astype(np.float64)

            constantError2 = 7.1 * 7.1
            self.sigmasq = self.imgarray + constantError2
            self.weight = 1.0/self.sigmasq

            # close file
            hdulist.close()

        else:
            self.inputHeader = {}
            extname = 'None'
            self.imgarray = inputImageArray.astype(np.float64)
            self.weight = 1.0/np.sqrt(self.imgarray)

        # setup starting Zernike array - start with the inputZernikeDict keyed by extension name
        # if wavefrontMap exists, then use that to fill all Zernike coefficients from zern5 to nZernikeSize (zern5 is iZ=3)

        # take this from the inputZernikeDict keyed by the value of extname in the header
        inputZernikeArray = self.fitDict["inputZernikeDict"][extname]

        # fill startingZernikeArray from the inputZernikeArray values
        startingZernikeArray = np.zeros(self.gFitFunc.nZernikeSize)
        for iZ in range(len(inputZernikeArray)):
            # starting from 1st used Zernike term
            startingZernikeArray[iZ] = inputZernikeArray[iZ]

        # if the wavefrontMap is present, overwrite the zernike terms after Focus (starting from 5 and up to nZernikeTerms) to values from Zemax built map
        if self.wavefrontMap != None:
            anotherZernikeArray = self.wavefrontMap.get(xDECam,yDECam,nZernikeFirst=5,nZernikeLast=self.paramDict["nZernikeTerms"])
            for iZ in range(5,self.paramDict["nZernikeTerms"]+1):
                startingZernikeArray[iZ-2] = anotherZernikeArray[iZ-5]

        self.startingParam[self.gFitFunc.ipar_ZernikeFirst:self.gFitFunc.ipar_ZernikeLast+1] = startingZernikeArray
        startingZernikeError = 0.1 * np.ones(self.gFitFunc.nZernikeSize)
        self.errorParam[self.gFitFunc.ipar_ZernikeFirst:self.gFitFunc.ipar_ZernikeLast+1] = startingZernikeError

        # set sign of Focus term
        # 3/25/2014 make the maximum Zern4 a setable parameter (for big donuts!)
        if np.sign(startingZernikeArray[2])==1:
            self.loParam[self.gFitFunc.ipar_ZernikeFirst+2] = 0.0
            self.hiParam[self.gFitFunc.ipar_ZernikeFirst+2] = self.fitDict["inputZ4max"]
        else:
            self.loParam[self.gFitFunc.ipar_ZernikeFirst+2] = -1 * self.fitDict["inputZ4max"]
            self.hiParam[self.gFitFunc.ipar_ZernikeFirst+2] = 0.0

        # set range for x,y tilts, to prevent runaway values
        # offset = 4 lambda F * zern_2 or zern_3 = 8.4 micron/ unit of z2,3
        # limit to +- 15 pixels = +-27 units, make it +-30
        self.loParam[self.gFitFunc.ipar_ZernikeFirst+0] = -30.0
        self.hiParam[self.gFitFunc.ipar_ZernikeFirst+0] = 30.0
        self.loParam[self.gFitFunc.ipar_ZernikeFirst+1] = -30.0
        self.hiParam[self.gFitFunc.ipar_ZernikeFirst+1] = 30.0

        # calculate staring values for nele and bkgd or take from paramDict
        background,nelectrons = calcStarting(self.imgarray,printLevel=self.paramDict["printLevel"])
        inputnEle = self.fitDict["inputnEle"]
        if inputnEle!=None:
            self.startingParam[self.gFitFunc.ipar_nEle] = inputnEle
            self.errorParam[self.gFitFunc.ipar_nEle] = np.sqrt(inputnEle)
        else:
            self.startingParam[self.gFitFunc.ipar_nEle] = nelectrons
            self.errorParam[self.gFitFunc.ipar_nEle] = np.sqrt(nelectrons)
        self.loParam[self.gFitFunc.ipar_nEle] = 0.2*self.startingParam[self.gFitFunc.ipar_nEle]
        self.hiParam[self.gFitFunc.ipar_nEle] = 4.0*self.startingParam[self.gFitFunc.ipar_nEle]

        inputbkgd = self.fitDict["inputbkgd"]
        if inputbkgd!=None:
            self.startingParam[self.gFitFunc.ipar_bkgd] = inputbkgd
            self.errorParam[self.gFitFunc.ipar_bkgd] = 1.0
            self.loParam[self.gFitFunc.ipar_bkgd] = 0.0
            self.hiParam[self.gFitFunc.ipar_bkgd] = inputbkgd*10.0
        else:
            self.startingParam[self.gFitFunc.ipar_bkgd] = background
            self.errorParam[self.gFitFunc.ipar_bkgd] = 1.0
            self.loParam[self.gFitFunc.ipar_bkgd] = 0.0
            self.hiParam[self.gFitFunc.ipar_bkgd] = background*10.0

        # rzero parameter, get from argument first then header or otherwise set to 0.2
            # March 25,2014 increase limit on rzero from 0.25 to 0.30
            # temporary - increase rzero limit to 0.50
        inputrzero = self.fitDict["inputrzero"]
        if inputrzero == None:
            if list(self.inputHeader.keys()).count("RZEROIN")>0:
                inputrzero = self.inputHeader["RZEROIN"]
            else:
                inputrzero = 0.125

        self.startingParam[self.gFitFunc.ipar_rzero] = inputrzero
        self.errorParam[self.gFitFunc.ipar_rzero] = 0.01
        self.loParam[self.gFitFunc.ipar_rzero] = 0.05
        self.hiParam[self.gFitFunc.ipar_rzero] = 0.50

        # Set starting values and step sizes for parameters
        # (note that one can redefine the parameters, so this method can be called multiple times)
        for ipar,aname in enumerate(self.gFitFunc.parNames):
            self.gMinuit.values[aname] = self.startingParam[ipar]
            if self.loParam[ipar]!=0 or self.hiParam[ipar]!=0 :
                self.gMinuit.limits[aname] = (self.loParam[ipar],self.hiParam[ipar])
            self.gMinuit.errors[aname] = self.errorParam[ipar]
            if self.paramStatusArray[ipar]==0:
                self.gMinuit.fixed[aname] = False
            else:
                self.gMinuit.fixed[aname] = True

        # do the Fit, and repeat as desired with different parameters fixed
        postfix = {0:"first",1:"second",2:"third",3:"fourth"}
        for iFit in range(self.paramDict["nFits"]):

            # fix parameters as desired
            fixParamArray = self.paramDict["fixedParamArray"+str(iFit+1)]
            self.updateFit(fixParamArray,xDECam,yDECam)
            self.doFit(self.fitDict)
            outputDict = self.outFit(postfix[iFit])

#        # if desired do the Wavefront fit next
#        if self.paramDict["doGridFit"]:
#
#            self.initWavefrontGrid(self.paramDict['spacing'])
#            self.initWavefrontFit()
#
#            #self.debugWavefrontFit()
#            #self.chisqvsparWavefrontFit(iPar=47)
#            self.calcderivWavefrontFit()
#
#            self.doWavefrontFit()
#            self.outWavefrontFit("gridfit")

        # return last output Dictionary
        return outputDict




    def updateFit(self,fixParamArray,xDECam,yDECam):
        # fix parameters as desired
        for ipar,aname in enumerate(self.gFitFunc.parNames):
            if fixParamArray[ipar] != self.paramStatusArray[ipar]:
                if fixParamArray[ipar]==0:
                    self.gMinuit.fixed[aname] = False
                elif fixParamArray[ipar]==1:
                    self.gMinuit.fixed[aname] = True
            self.paramStatusArray[ipar] = fixParamArray[ipar]

        # set x,y DECam values
        self.gFitFunc.setXYDECam(xDECam,yDECam)

        # reset counters
        self.gFitFunc.nCallsCalcAll = 0
        self.gFitFunc.nCallsCalcDerivative = 0




    def doFit(self,kwargs):

        # print out
        print(tab.tabulate(*self.gMinuit.params.to_table()))

        # set max iterations here

        # start timer
        self.startingtime = time.time()

        # do the fit
        self.gMinuit.migrad()

        # add option to not calculate the gradient...
#        if not kwargs['skipGradient']:
#            # tell Minuit we have derivatives, don't check anymore as long as rzero is fixed
#            if self.paramStatusArray[self.gFitFunc.ipar_rzero]==1 :
#                self.gFitFunc.setCalcRzeroDerivativeFalse()
#                arglist[0] = 1  # =1 means never check the gradient
#                self.gMinuit.mnexcm( "SET GRADIENT", arglist, 1, ierflg )
#            else:
#                self.gFitFunc.setCalcRzeroDerivativeTrue()
#                arglist[0] = 1  # =1 means never check the gradient (March 12, 2015, rzero derivative code implemented)
                #            arglist[0] = 0  # =0 means to check gradient each time
#                self.gMinuit.mnexcm( "SET GRADIENT", arglist, 1, ierflg )
#        else:
#            self.gFitFunc.setCalcRzeroDerivativeFalse()


        # done, check elapsed time
        firsttime = time.time()
        self.deltatime = firsttime - self.startingtime
        if self.paramDict["printLevel"]>=1:
            print('donutfit: Elapsed time fit = ',self.deltatime)

        # number of calls
        if self.paramDict["printLevel"]>=1:
            print('donutfit: Number of CalcAll calls = ',self.gFitFunc.nCallsCalcAll)
            print('donutfit: Number of CalcDerivative calls = ',self.gFitFunc.nCallsCalcDerivative)

        # print results
        print(tab.tabulate(*self.gMinuit.params.to_table()))


    def outFit(self,postfix,identifier=""):

        # get fit details from iminuit
        amin = self.gMinuit.fmin.fval
        edm = self.gMinuit.fmin.edm
        errdef = self.gMinuit.fmin.errordef
        nvpar = self.gMinuit.nfit
        nparx = self.gMinuit.npar
        icstat = int(self.gMinuit.fmin.is_valid) + 2*int(self.gMinuit.fmin.has_accurate_covar)
        dof = pow(self.gFitFunc._nPixels,2) - nvpar

        if self.paramDict["printLevel"]>=1:
            mytxt = "amin = %.3f, edm = %.3f,   effdef = %.3f,   nvpar = %d,  nparx = %d, icstat = %d " % (amin,edm,errdef,nvpar,nparx,icstat)
            print('donutfit: ',mytxt)

        # get fit values and errors
        self.paramArray = self.gMinuit.values
        self.paramErrArray = self.gMinuit.errors

        # printout Zernike parameters in a convenient format
        if self.paramDict["printLevel"]>=2:
            print(""" "[ """, end=' ')
            for ipar in range(self.gFitFunc.ipar_ZernikeFirst,self.gFitFunc.ipar_ZernikeLast):
                print(self.paramArray[ipar],",", end=' ')
            print(self.paramArray[self.gFitFunc.ipar_ZernikeLast],""" ]" """)

        #copy input header information from input file here except for Standard header stuff
        stdHeaderDict= {'SIMPLE':0,'BITPIX':0,'NAXIS':0,'NAXIS1':0,'NAXIS2':0,'EXTEND':0}
        outputHeaderDict = OrderedDict()

        for key in list(self.inputHeader.keys()):
            if not list(stdHeaderDict.keys()).count(key)>0:
                outputHeaderDict[key] = self.inputHeader[key]

        # fill output Dictionary
        outputDict = OrderedDict()

        outputDict["CHI2"] = amin
        outputDict["DOF"] = dof
        outputDict["FITSTAT"] = icstat
        outputDict["CLKTIME"] = self.deltatime
        outputDict["NCALCALL"] = self.gFitFunc.nCallsCalcAll
        outputDict["NCALCDER"] = self.gFitFunc.nCallsCalcDerivative

        for ipar in range(self.gFitFunc.npar):
            outputDict[self.gFitFunc.parNames[ipar]] = float(self.paramArray[ipar])
        for ipar in range(self.gFitFunc.npar):
            outputDict[self.gFitFunc.parNames[ipar]+"E"] = float(self.paramErrArray[ipar])

        # make a single output file, with multiple extensions
        # Extension 1:  Calculated Image
        # Extension 2:  Original Image
        # Extension 3:  Difference  (if desired)
        # Extension 4:  Chi2        (if desired)
        # Extension 5:  Wavefront   (if desired)

        hduListOutput = pyfits.HDUList()
        primaryOutput = pyfits.PrimaryHDU()
        primaryHeader = primaryOutput.header

        # fill primary header both with input Header and fit results
        for key in outputHeaderDict:
            primaryHeader[key] = outputHeaderDict[key]
        for key in outputDict:
            primaryHeader[key] = outputDict[key]
        hduListOutput.append(primaryOutput)

        # calculated Donut
        calcHdu = pyfits.ImageHDU(self.gFitFunc.getvImage(),name="Model_Image")
        calcHeader = calcHdu.header
        for key in outputHeaderDict:
            calcHeader[key] = outputHeaderDict[key]
        for key in outputDict:
            calcHeader[key] = outputDict[key]
        hduListOutput.append(calcHdu)

        # original image
        imageHdu = pyfits.ImageHDU(self.imgarray,name="Original_Image")
        imageHeader = imageHdu.header
        for key in outputHeaderDict:
            imageHeader[key] = outputHeaderDict[key]
        hduListOutput.append(imageHdu)

        # diff Donut - Calc
        if self.paramDict["outputDiff"]:
            diffHdu = pyfits.ImageHDU(self.imgarray-self.gFitFunc.getvImage(),name="Image-Model")
            diffHeader = diffHdu.header
            for key in outputHeaderDict:
                imageHeader[key] = outputHeaderDict[key]
            hduListOutput.append(diffHdu)

        # Chi2 Donut-Calc
        if self.paramDict["outputChi2"]:
            chi2Hdu = pyfits.ImageHDU(self.pullsq,name="Chi2 Image-Model")
            chi2Header = chi2Hdu.header
            for key in outputHeaderDict:
                imageHeader[key] = outputHeaderDict[key]
            hduListOutput.append(chi2Hdu)


        # Zernike Wavefront map
        if self.paramDict["outputWavefront"]:
            waveHdu = pyfits.ImageHDU(self.gFitFunc.getvPupilMask()*self.gFitFunc.getvPupilWaveZernike(),name="Total_Wavefront")
            waveHeader = waveHdu.header
            hduListOutput.append(waveHdu)

        # file names for output are
        # outputPrefix + identifier
        if identifier!="" :
            outName = self.fitDict["outputPrefix"] + "." + identifier + "." + postfix
        else:
            outName = self.fitDict["outputPrefix"] + "." + postfix

        # write out fits file
        outFile =  outName + ".donut.fits"
        if self.paramDict["printLevel"]>=1:
            hduListOutput.info()
        hduListOutput.writeto(outFile,overwrite=True)

        # add info from input Header for return
        outputDict.update(outputHeaderDict)
        resultDict = {}
        resultDict.update(outputDict)
        return resultDict


    def initWavefrontGrid(self,spacing):

        self.gFitFunc.initWavefrontGrid(spacing)

    def initWavefrontFit(self):

        # make a new gMinit object
        nGrid = self.gFitFunc.getnGrid()
        self.wMinuit = ROOT.TMinuit(nGrid)
        self.wMinuit.SetFCN(self.wavechisq)

        # initialize array of paramters from last iteration
        self.lastParArr = np.zeros(nGrid)

        # arglist is for the parameters in Minuit commands
        arglist = array( 'd', 10*[0.] )
        ierflg =ctypes.c_int(1982)#L ROOT.Long(1982)

        # set the definition of 1sigma
        arglist[0] = 1.0
        self.wMinuit.mnexcm( "SET ERR", arglist, 1, ierflg )

        # turn off Warnings
        arglist[0] = 0
        self.wMinuit.mnexcm("SET NOWARNINGS", arglist,0,ierflg)

        # set printlevel
        arglist[0] = self.paramDict["printLevel"]
        self.wMinuit.mnexcm("SET PRINTOUT", arglist,1,ierflg)

        # do initial setup of Minuit parameters
        startingwParam = np.zeros(nGrid)
        errorwParam = 0.01 * np.ones(nGrid)
        maxwavevalue = 0.5
        lowParam = -1.0 * np.ones(nGrid) * maxwavevalue
        hiwParam = np.ones(nGrid) * maxwavevalue
        wparamStatusArray = np.zeros(nGrid)   # store =0 Floating, =1 Fixed

        # Set starting values and step sizes for parameters
        wparamNames = []
        for ipar in range(nGrid):
            wparamNames.append("Grid_%d" % (ipar))
            self.wMinuit.DefineParameter(ipar,"Grid_%d" % (ipar), startingwParam[ipar], errorwParam[ipar],lowParam[ipar],hiwParam[ipar])

        self.nCallsWavefit = 0


    def wavechisq(self,npar,gin,f,par,iflag):

        # convert par to a np array
        nGrid = self.gFitFunc.getnGrid()
        parArr = np.zeros(nGrid)  # convert npar to an integer...
        for ipar in range(nGrid):
            parArr[ipar] = par[ipar]

        # call donutengine to calculate image
        self.gFitFunc.nCallsCalcAll +=1
        self.gFitFunc.makeWavefrontGrid(parArr)  # fill deltaWFM from Minuit parameters
        self.gFitFunc.calcAll(self.gFitFunc.getParCurrent())   # calculate Image  # still a bug?? ### BUG - needs fixing to use most recent regular parameters from regular fit

        # write out the image here and quit...
        if False:
            vImage = self.gFitFunc.getvImage()
            ftemp = open('image_%d.npy' % self.gFitFunc.nCallsCalcAll, 'wb')
            np.save(ftemp,vImage)

        # compare to calculated image
        diff = self.imgarray - self.gFitFunc.getvImage()
        self.pullsq = diff*diff/self.sigmasq
        chisquared = self.pullsq.sum()

        # how many parameters changed from last time...
        nChange = 0
        parChange = []
        for ipar in range(nGrid):
            if parArr[ipar] != self.lastParArr[ipar]:
                nChange = nChange+1
                parChange.append([ipar,parArr[ipar],self.lastParArr[ipar]])
            self.lastParArr[ipar] = parArr[ipar]

        # printout
        #if self.paramDict["printLevel"]>=2:
        print('donutfit: wavechisq Chi2 = ',chisquared, ' Nchanged Par ',nChange,' ',parChange)
        #else:
        #    if nChange>2 :
        #        print('donutfit: wavechisq Chi2 = ',chisquared, ' Nchanged Par ',nChange)

        # return result
        f[0] = chisquared

        # Derivative calculation!
#        if iflag==2 :

            # not currently called in calcAll
#            self.gFitFunc.calcGridDerivatives(self.imgarray,self.weight)
#            dChi2dgrid = self.gFitFunc.getGridDerivatives()

#            gin.SetSize(nGrid)  # need to handle root bug
            #
            # fill gin with Derivatives
            #
#            print("CalcGridDerivatives:")
#            for i in range(nGrid):
#                gin[i] = dChi2dgrid[i]
#                #print(i,gin[i])

    def chisqvsparWavefrontFit(self,iPar):

        # plot chi2 vs. parameter
        parArr = np.zeros(self.gFitFunc.getnGrid())

        parvals = np.arange(-0.5,0.5,0.01)
        chivals = []

        for val in parvals:
            parArr[iPar] = val
            self.gFitFunc.makeWavefrontGrid(parArr)
            self.gFitFunc.calcAll(self.gFitFunc.getParCurrent())   # calculate Image

            # compare to calculated image
            diff = self.imgarray - self.gFitFunc.getvImage()
            self.pullsq = diff*diff/self.sigmasq
            chivals.append(self.pullsq.sum())

        print(chivals)

    def calcderivWavefrontFit(self):

        # plot derivative vs. Grid parameter
        nGrid = self.gFitFunc.getnGrid()
        parArr = np.zeros(nGrid)
        deriv = np.zeros(nGrid)
        delta = 0.05

        iiGrid = self.gFitFunc.getviiGrid()
        jjGrid = self.gFitFunc.getvjjGrid()
        nbinGrid = self.gFitFunc._nbinGrid
        deriv2d = np.zeros((nbinGrid,nbinGrid))

        self.gFitFunc.makeWavefrontGrid(parArr)
        self.gFitFunc.calcAll(self.gFitFunc.getParCurrent())   # calculate Image

        # compare to calculated image
        diff = self.imgarray - self.gFitFunc.getvImage()
        pullsq = diff*diff/self.sigmasq
        chisq0 = pullsq.sum()

        for iGrid in range(nGrid):
            parArr[iGrid] = delta
            self.gFitFunc.makeWavefrontGrid(parArr)
            self.gFitFunc.calcAll(self.gFitFunc.getParCurrent())   # calculate Image

            # compare to calculated image
            diff = self.imgarray - self.gFitFunc.getvImage()
            pullsq = diff*diff/self.sigmasq
            derivative = (pullsq.sum() - chisq0)/delta
            deriv[iGrid] = derivative
            deriv2d[jjGrid[iGrid],iiGrid[iGrid]] = deriv[iGrid]

            # if log(abs(deriv))< 2.0 then fix Paramter
            if np.log10(np.abs(derivative))<2.0 :
                self.wMinuit.FixParameter(iGrid)

            # reset
            parArr[iGrid] = 0.0


        hduListOutput = pyfits.HDUList()
        hduListOutput.append(pyfits.ImageHDU(deriv2d))
        hduListOutput.writeto("grid-deriv.fits",overwrite=True)



    def debugWavefrontFit(self):

        hduListOutput = pyfits.HDUList()

        # loop over parameters - change one at a time, make the fineWavefrontGrid and plot both grids
        for i in range(self.gFitFunc.getnGrid()):
            parArr = np.zeros(self.gFitFunc.getnGrid())
            parArr[i] = 0.1
            self.gFitFunc.makeWavefrontGrid(parArr)  # fill deltaWFM from Minuit parameters

            # Wavefront Grid map
            hduListOutput.append(pyfits.ImageHDU(self.gFitFunc.getvDeltaWFM().copy()))

            # Coarse Wavefront Grid map
            hduListOutput.append(pyfits.ImageHDU(self.gFitFunc.getvCoarseWFM().copy()))

        hduListOutput.writeto("/u/ec/roodman/kipacdisk/grid-debug.fits",overwrite=True)


    def doWavefrontFit(self):

        # arglist is for the parameters in Minuit commands
        arglist = array( 'd', 10*[0.] )
        ierflg = ctypes.c_int(1982) #L ROOT.Long

        # tell Minuit we have derivatives, don't check anymore as long as rzero is fixed
        #arglist[0] = 1  # =1 means never check the gradient
        arglist[0] = 0  # =0 means to check gradient each time
        self.wMinuit.mnexcm( "SET GRADIENT", arglist, 1, ierflg )

        # tell Minuit to use strategy for fastest fits
        arglist[0] = 0  # was 1
        self.wMinuit.mnexcm( "SET STRATEGY", arglist, 1, ierflg )

        # start timer
        self.startingtime = time.time()

        # Now ready for minimization step
        arglist[0] = self.paramDict["maxIterations"]
        arglist[1] = 1000.    # tolerance, default is 0.1
        self.wMinuit.mnexcm( "MIGRAD", arglist, 2, ierflg )

        # done, check elapsed time
        firsttime = time.time()
        self.deltatime = firsttime - self.startingtime
        if self.paramDict["printLevel"]>=1:
            print('donutfit: Elapsed time fit = ',self.deltatime)

        # number of calls
        if self.paramDict["printLevel"]>=1:
            print('donutfit wavefront: Number of CalcAll calls = ',self.gFitFunc.nCallsCalcAll)
            print('donutfit wavefront: Number of CalcDerivative calls = ',self.gFitFunc.nCallsCalcDerivative)

    def outWavefrontFit(self,postfix,identifier=""):

        # get more fit details from MINUIT
        amin, edm, errdef = ctypes.c_double(0.18), ctypes.c_double(0.19), ctypes.c_double(0.20)
        nvpar, nparx, icstat = ctypes.c_int(1983), ctypes.c_int(1984), ctypes.c_int(1985)
        self.wMinuit.mnstat( amin, edm, errdef, nvpar, nparx, icstat )
        dof = pow(self.gFitFunc._nPixels,2) - nvpar.value
        if self.paramDict["printLevel"]>=1:
            mytxt = "amin = %.3f, edm = %.3f,   effdef = %.3f,   nvpar = %.3f,  nparx = %.3f, icstat = %.3f " % (amin.value,edm.value,errdef.value,nvpar.value,nparx.value,icstat.value)
            print('donutfit wavefront: ',mytxt)

        # get fit values and errors
        aVal = ctypes.c_double(0.21)
        errVal = ctypes.c_double(0.22)
        self.paramArray = np.zeros(self.gFitFunc.getnGrid())
        self.paramErrArray = np.zeros(self.gFitFunc.getnGrid())
        nGrid = self.gFitFunc.getnGrid()
        for ipar in range(nGrid):
            self.wMinuit.GetParameter(ipar,aVal,errVal)
            self.paramArray[ipar] = aVal.value
            if errVal.value < 1e9 :
                self.paramErrArray[ipar] = errVal.value
            else:
                self.paramErrArray[ipar] = 0.0

        #copy input header information from input file here except for Standard header stuff
        stdHeaderDict= {'SIMPLE':0,'BITPIX':0,'NAXIS':0,'NAXIS1':0,'NAXIS2':0,'EXTEND':0}
        outputHeaderDict = OrderedDict()

        for key in list(self.inputHeader.keys()):
            if not list(stdHeaderDict.keys()).count(key)>0:
                outputHeaderDict[key] = self.inputHeader[key]

        # fill output Dictionary
        outputDict = OrderedDict()

        outputDict["CHI2"] = float(amin.value)
        outputDict["DOF"] = dof
        outputDict["FITSTAT"] = int(icstat.value)
        outputDict["CLKTIME"] = self.deltatime
        outputDict["NCALCALL"] = self.gFitFunc.nCallsCalcAll
        outputDict["NCALCDER"] = self.gFitFunc.nCallsCalcDerivative

        # do I use this anywhere?
        #for ipar in range(nGrid):
        #    outputDict[self.gFitFunc.parNames[ipar]] = float(self.paramArray[ipar])
        #for ipar in range(nGrid):
        #    outputDict[self.gFitFunc.parNames[ipar]+"E"] = float(self.paramErrArray[ipar])

        # make a single output file, with multiple extensions
        # Extension 1:  Calculated Image
        # Extension 2:  Original Image
        # Extension 3:  Difference  (if desired)
        # Extension 4:  Chi2        (if desired)
        # Extension 5:  Wavefront   (if desired)
        # Extension 6:  Fine Wavefront Grid
        # Extension 7:  Coarse Wavefront Grid

        hduListOutput = pyfits.HDUList()
        primaryOutput = pyfits.PrimaryHDU()
        primaryHeader = primaryOutput.header

        # fill primary header both with input Header and fit results
        for key in outputHeaderDict:
            primaryHeader[key] = outputHeaderDict[key]
        for key in outputDict:
            primaryHeader[key] = outputDict[key]
        hduListOutput.append(primaryOutput)

        # calculated Donut
        calcHdu = pyfits.ImageHDU(self.gFitFunc.getvImage())
        calcHeader = calcHdu.header
        for key in outputHeaderDict:
            calcHeader[key] = outputHeaderDict[key]
        for key in outputDict:
            calcHeader[key] = outputDict[key]
        hduListOutput.append(calcHdu)

        # original image
        imageHdu = pyfits.ImageHDU(self.imgarray)
        imageHeader = imageHdu.header
        for key in outputHeaderDict:
            imageHeader[key] = outputHeaderDict[key]
        hduListOutput.append(imageHdu)

        # diff Donut - Calc
        if self.paramDict["outputDiff"]:
            diffHdu = pyfits.ImageHDU(self.imgarray-self.gFitFunc.getvImage())
            diffHeader = diffHdu.header
            for key in outputHeaderDict:
                imageHeader[key] = outputHeaderDict[key]
            hduListOutput.append(diffHdu)

        # Chi2 Donut-Calc
        if self.paramDict["outputChi2"]:
            chi2Hdu = pyfits.ImageHDU(self.pullsq)
            chi2Header = chi2Hdu.header
            for key in outputHeaderDict:
                imageHeader[key] = outputHeaderDict[key]
            hduListOutput.append(chi2Hdu)


        # Wavefront map
        if self.paramDict["outputWavefront"]:
            waveHdu = pyfits.ImageHDU(self.gFitFunc.getvPupilMask()*self.gFitFunc.getvPupilWaveZernike())
            waveHeader = waveHdu.header
            hduListOutput.append(waveHdu)

        # Wavefront Grid map
        gridHdu = pyfits.ImageHDU(self.gFitFunc.getvDeltaWFM())
        gridHeader = gridHdu.header
        hduListOutput.append(gridHdu)

        # Coarse Wavefront Grid map
        gridHdu = pyfits.ImageHDU(self.gFitFunc.getvCoarseWFM())
        gridHeader = gridHdu.header
        hduListOutput.append(gridHdu)

        # file names for output are
        # outputPrefix + identifier
        if identifier!="" :
            outName = self.fitDict["outputPrefix"] + "." + identifier + "." + postfix
        else:
            outName = self.fitDict["outputPrefix"] + "." + postfix

        # write out fits file
        outFile =  outName + ".donut.fits"
        if self.paramDict["printLevel"]>=1:
            hduListOutput.info()
        hduListOutput.writeto(outFile,overwrite=True)

        # add info from input Header for return
        outputDict.update(outputHeaderDict)
        resultDict = {}
        resultDict.update(outputDict)
        return resultDict
