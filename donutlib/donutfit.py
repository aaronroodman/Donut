import time
import os
import numpy
import sys
# OrderedDict is a feature of 2.7 and beyond only
from collections import OrderedDict
from astropy.io import fits as pyfits
from array import array
import pdb

import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True
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
                          "calcRzeroDerivative":True}

        # search for key in inputDict, change defaults
        self.paramDict.update(inputDict)

        # setup the fit engine
        self.gFitFunc = donutengine(**self.paramDict)
        
        # need dummy versions before calling self.chisq
        self.imgarray = numpy.zeros(1)
        self.weight = numpy.ones(1)
        self.sigmasq = numpy.ones(1)

        # get decam info
        self.decamInfo = decaminfo()

        # setup MINUIT
        self.gMinuit = ROOT.TMinuit(self.gFitFunc.npar)
        self.gMinuit.SetFCN( self.chisq )

        # arglist is for the parameters in Minuit commands
        arglist = array( 'd', 10*[0.] )
        ierflg = ROOT.Long(1982)

        # set the definition of 1sigma 
        arglist[0] = 1.0
        self.gMinuit.mnexcm( "SET ERR", arglist, 1, ierflg )

        # turn off Warnings
        arglist[0] = 0
        self.gMinuit.mnexcm("SET NOWARNINGS", arglist,0,ierflg)

        # set printlevel
        arglist[0] = self.paramDict["printLevel"]
        self.gMinuit.mnexcm("SET PRINTOUT", arglist,1,ierflg)

        # do initial setup of Minuit parameters

        # status/limit arrays for Minuit parameters
        self.startingParam = numpy.zeros(self.gFitFunc.npar)
        self.errorParam = numpy.ones(self.gFitFunc.npar)
        self.loParam = numpy.zeros(self.gFitFunc.npar)
        self.hiParam = numpy.zeros(self.gFitFunc.npar)
        self.paramStatusArray = numpy.zeros(self.gFitFunc.npar)   # store =0 Floating, =1 Fixed

        # Set starting values and step sizes for parameters
        # (note that one can redefine the parameters, so this method can be called multiple times)
        for ipar in range(self.gFitFunc.npar):
            self.gMinuit.DefineParameter(ipar,self.gFitFunc.parNames[ipar],self.startingParam[ipar],self.errorParam[ipar],self.loParam[ipar],self.hiParam[ipar])


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
                        "inputZ4max":100.0}   #Sept 22, 2014 change default back to 100

        # update Dictionary with inputs
        self.fitDict.update(inputFitDict)

        # if desired print dict
        if self.paramDict["printLevel"]>=2:
            print self.fitDict
        
        # load image either from a file or from an input array
        if self.fitDict["inputImageArray"] ==None:
            
            # check for fits in file name
            fileName = self.fitDict["inputFile"]
            if fileName[-5:] != ".fits" :
                fileName = fileName + ".fits"

            # get the input file header
            hdulist = pyfits.open(fileName)
            self.inputHeader = hdulist[0].header

            # get the extension name from the header
            if self.inputHeader.keys().count("EXTNAME")>0:
                extname = self.inputHeader["EXTNAME"]
                if extname == "":
                    extname = 'None'
            else:
                extname = 'None'

            # also get the IX,IY values from the header
            if self.inputHeader.keys().count("IX")>0:
                ix = self.inputHeader["IX"]
            else:
                ix = 0.
                
            if self.inputHeader.keys().count("IY")>0:
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
            self.imgarray = hdulist[0].data.astype(numpy.float64)
                
            constantError2 = 7.1 * 7.1 
            self.sigmasq = self.imgarray + constantError2
            self.weight = 1.0/self.sigmasq

            # close file
            hdulist.close()

        else:
            self.inputHeader = {}
            extname = 'None'
            self.imgarray = inputImageArray.astype(numpy.float64)
            self.weight = 1.0/numpy.sqrt(self.imgarray)
       
        # setup starting Zernike array
        # take this from the inputZernikeDict keyed by the value of extname in the header
        inputZernikeArray = self.fitDict["inputZernikeDict"][extname]
        startingZernikeArray = numpy.zeros(self.gFitFunc.nZernikeSize)
        for iZ in range(len(inputZernikeArray)):
            # starting from 1st used Zernike term
            startingZernikeArray[iZ] = inputZernikeArray[iZ]

        self.startingParam[self.gFitFunc.ipar_ZernikeFirst:self.gFitFunc.ipar_ZernikeLast+1] = startingZernikeArray
        startingZernikeError = 0.1 * numpy.ones(self.gFitFunc.nZernikeSize)
        self.errorParam[self.gFitFunc.ipar_ZernikeFirst:self.gFitFunc.ipar_ZernikeLast+1] = startingZernikeError

        # set sign of Focus term 
        # 3/25/2014 make the maximum Zern4 a setable parameter (for big donuts!)
        if numpy.sign(startingZernikeArray[2])==1:
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
            self.errorParam[self.gFitFunc.ipar_nEle] = numpy.sqrt(inputnEle)
        else:
            self.startingParam[self.gFitFunc.ipar_nEle] = nelectrons
            self.errorParam[self.gFitFunc.ipar_nEle] = numpy.sqrt(nelectrons)
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
            if self.inputHeader.keys().count("RZEROIN")>0:
                inputrzero = self.inputHeader["RZEROIN"]
            else:
                inputrzero = 0.125
            
        self.startingParam[self.gFitFunc.ipar_rzero] = inputrzero
        self.errorParam[self.gFitFunc.ipar_rzero] = 0.01
        self.loParam[self.gFitFunc.ipar_rzero] = 0.05
        self.hiParam[self.gFitFunc.ipar_rzero] = 0.50

        # Set starting values and step sizes for parameters
        # (note that one can redefine the parameters, so this method can be called multiple times)
        for ipar in range(self.gFitFunc.npar):
            self.gMinuit.DefineParameter(ipar,self.gFitFunc.parNames[ipar],self.startingParam[ipar],self.errorParam[ipar],self.loParam[ipar],self.hiParam[ipar])

        # do the Fit, and repeat as desired with different parameters fixed
        postfix = {0:"first",1:"second",2:"third",3:"fourth"}
        for iFit in range(self.paramDict["nFits"]):

            # fix parameters as desired
            fixParamArray = self.paramDict["fixedParamArray"+str(iFit+1)]
            self.updateFit(fixParamArray,xDECam,yDECam)
            self.doFit()
            outputDict = self.outFit(postfix[iFit])

        # return last output Dictionary
        return outputDict

    def updateFit(self,fixParamArray,xDECam,yDECam):
        # fix parameters as desired
        for ipar in range(self.gFitFunc.npar):
            if fixParamArray[ipar] != self.paramStatusArray[ipar]:
                if fixParamArray[ipar]==0:
                    self.gMinuit.Release(ipar)
                elif fixParamArray[ipar]==1:
                    self.gMinuit.FixParameter(ipar)                    
            self.paramStatusArray[ipar] = fixParamArray[ipar]

        # set x,y DECam values
        self.gFitFunc.setXYDECam(xDECam,yDECam)

        # reset counters
        self.gFitFunc.nCallsCalcAll = 0
        self.gFitFunc.nCallsCalcDerivative = 0


    def chisq(self,npar, gin, f, par, iflag ):

        # convert par to a numpy array
        parArr = numpy.zeros(self.gFitFunc.npar)
        for ipar in range(self.gFitFunc.npar):
            parArr[ipar] = par[ipar]
        
        # call donutengine to calculate image
        self.gFitFunc.calcAll(parArr)

        # compare to calculated image
        diff = self.imgarray - self.gFitFunc.getvImage()
        self.pullsq = diff*diff/self.sigmasq   
        chisquared = self.pullsq.sum()

        # printout
        if self.paramDict["printLevel"]>=2:
            print 'donutfit: Chi2 = ',chisquared

        # save parameters for next iteration
        self.gFitFunc.savePar()

        # return result    
        f[0] = chisquared

        if iflag==2 :

            # not currently called in calcAll
            self.gFitFunc.calcDerivatives(self.imgarray,self.weight)
            dChi2dpar = self.gFitFunc.getDerivatives()
            
            gin.SetSize(self.gFitFunc.npar)  # need to handle root bug
            #
            # fill gin with Derivatives
            #            
            for i in range(self.gFitFunc.npar):
                gin[i] = dChi2dpar[i]

        #print "donutfit.chisq: ",iflag,f[0]
        #print "                ",parArr
        #if iflag==2:
        #    print "                ",dChi2dpar

    def doFit(self):
        
        # arglist is for the parameters in Minuit commands
        arglist = array( 'd', 10*[0.] )
        ierflg = ROOT.Long(1982)

        # tell Minuit we have derivatives, don't check anymore as long as rzero is fixed
        if self.paramStatusArray[self.gFitFunc.ipar_rzero]==1 :
            self.gFitFunc.setCalcRzeroDerivativeFalse()
            arglist[0] = 1  # =1 means never check the gradient
            self.gMinuit.mnexcm( "SET GRADIENT", arglist, 1, ierflg )
        else:
            self.gFitFunc.setCalcRzeroDerivativeTrue()
            arglist[0] = 1  # =1 means never check the gradient (March 12, 2015, rzero derivative code implemented)
#            arglist[0] = 0  # =0 means to check gradient each time
            self.gMinuit.mnexcm( "SET GRADIENT", arglist, 1, ierflg )

        # tell Minuit to use strategy for fastest fits
        arglist[0] = 0  # was 1
        self.gMinuit.mnexcm( "SET STRATEGY", arglist, 1, ierflg )
                
        # start timer
        self.startingtime = time.clock()

        # Now ready for minimization step
        self.gMinuit.SetMaxIterations(self.paramDict["maxIterations"])
        self.gMinuit.Migrad()

        # done, check elapsed time
        firsttime = time.clock()
        self.deltatime = firsttime - self.startingtime
        if self.paramDict["printLevel"]>=1:
            print 'donutfit: Elapsed time fit = ',self.deltatime

        # number of calls
        if self.paramDict["printLevel"]>=1:
            print 'donutfit: Number of CalcAll calls = ',self.gFitFunc.nCallsCalcAll
            print 'donutfit: Number of CalcDerivative calls = ',self.gFitFunc.nCallsCalcDerivative

    def outFit(self,postfix,identifier=""):

        # get more fit details from MINUIT
        amin, edm, errdef = ROOT.Double(0.18), ROOT.Double(0.19), ROOT.Double(0.20)
        nvpar, nparx, icstat = ROOT.Long(1983), ROOT.Long(1984), ROOT.Long(1985)
        self.gMinuit.mnstat( amin, edm, errdef, nvpar, nparx, icstat )
        dof = pow(self.gFitFunc._nPixels,2) - nvpar
        if self.paramDict["printLevel"]>=1:
            mytxt = "amin = %.3f, edm = %.3f,   effdef = %.3f,   nvpar = %.3f,  nparx = %.3f, icstat = %.3f " % (amin,edm,errdef,nvpar,nparx,icstat)
            print 'donutfit: ',mytxt

        # get fit values and errors
        aVal = ROOT.Double(0.21)
        errVal = ROOT.Double(0.22)
        self.paramArray = numpy.zeros(self.gFitFunc.npar)
        self.paramErrArray = numpy.zeros(self.gFitFunc.npar)
        for ipar in range(self.gFitFunc.npar):
            self.gMinuit.GetParameter(ipar,aVal,errVal)
            self.paramArray[ipar] = aVal
            if errVal < 1e9 :
                self.paramErrArray[ipar] = errVal
            else:
                self.paramErrArray[ipar] = 0.0

        # printout parameters in a convenient format
        if self.paramDict["printLevel"]>=1:
            print """ "[ """,
            for ipar in range(self.gFitFunc.ipar_ZernikeFirst,self.gFitFunc.ipar_ZernikeLast):
                print self.paramArray[ipar],",",
            print self.paramArray[self.gFitFunc.ipar_ZernikeLast],""" ]" """

        #copy input header information from input file here except for Standard header stuff
        stdHeaderDict= {'SIMPLE':0,'BITPIX':0,'NAXIS':0,'NAXIS1':0,'NAXIS2':0,'EXTEND':0}
        try:
            if sys.version_info.minor>=7:
                outputHeaderDict = OrderedDict()
        except:
            outputHeaderDict = {}

        for key in self.inputHeader.keys():
            if not stdHeaderDict.keys().count(key)>0:
                outputHeaderDict[key] = self.inputHeader[key]

        # fill output Dictionary
        try:
            if sys.version_info.minor>=7:
                outputDict = OrderedDict()
        except:
            outputDict = {}

        outputDict["CHI2"] = float(amin)
        outputDict["DOF"] = int(dof)
        outputDict["FITSTAT"] = int(icstat)
        outputDict["CLKTIME"] = self.deltatime
        outputDict["NCALCALL"] = self.gFitFunc.nCallsCalcAll
        outputDict["NCALCDER"] = self.gFitFunc.nCallsCalcDerivative
        outputDict["DOF"] = dof
        
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
        hduListOutput.writeto(outFile,clobber=True)

        # add info from input Header for return
        outputDict.update(outputHeaderDict)
        resultDict = {}
        resultDict.update(outputDict)
        return resultDict




