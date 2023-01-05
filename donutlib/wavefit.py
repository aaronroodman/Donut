import time
import numpy as np
import cv2
from collections import OrderedDict
from astropy.io import fits as pyfits
import tabulate as tab

from iminuit import Minuit
from donutlib.donutengine import donutengine

class wavefit(object):
    """ donutfit is a class used to fit the wavefront of donuts, using donutengine and MINUIT, for the DES experiment

    Aaron Roodman (C) SLAC National Accelerator Laboratory, Stanford University 2019.
    """

    def __init__(self,dfObject,**inputDict):
        # init contains all initializations which are done only once for all fits
        # parameters in fixParamArray1 are nEle,rzero,bkgd,Z2,Z3,Z4,....Z11

        self.paramDict = {"outputPrefix":"test","printLevel":2,"maxIterations":None,"spacing":16,"maxwavevalue":0.5,"tolerance":0.1}

        # search for key in inputDict, change defaults
        self.paramDict.update(inputDict)

        # save the input fit engine and donut engine
        self.df = dfObject
        self.gFitFunc = dfObject.gFitFunc

        # get needed info from Zernike fit
        self.parCurrent = self.gFitFunc.getParCurrent()
        bkg = self.parCurrent[self.gFitFunc.ipar_bkgd]
        self.pupilMask = self.gFitFunc.getvPupilMask()

        self._nbin = self.gFitFunc._nbin
        self._nPixels = self.gFitFunc._nPixels

        # points for coarse and fine grids on pupil
        self._spacing = self.paramDict["spacing"]
        self.coarsegrid = []
        self.finegrid = []
        for i,ii in enumerate(range(0,self._nbin,self._spacing)):
            for j,jj in enumerate(range(0,self._nbin,self._spacing)):
                if self.pupilMask[jj,ii] > 0:
                    self.coarsegrid.append([j,i])
                    self.finegrid.append([jj,ii])

        # Free parameters
        self.npar = len(self.finegrid)

        # get weights and image
        constantError2 = 7.1 * 7.1
        self.sigmasq = self.df.imgarray + constantError2
        self.weight = 1.0/self.sigmasq
        self.imgarrayc = self.df.imgarray

        # do initial setup of Minuit parameters
        self.setupCoarseGrid(values=np.zeros(self.npar))

        # setup MINUIT
        self.gMinuit = Minuit(self.chisq,self.startingParam,name=self.parNames)
        self.gMinuit.strategy = 0
        self.gMinuit.tol = self.paramDict["tolerance"]
        self.gMinuit.errordef = Minuit.LEAST_SQUARES

        # set parameter values,limits, errors
        for ipar,aname in enumerate(self.parNames):
            self.gMinuit.values[aname] = self.startingParam[ipar]
            if self.loParam[ipar]!=0 or self.hiParam[ipar]!=0 :
                self.gMinuit.limits[aname] = (self.loParam[ipar],self.hiParam[ipar])
            self.gMinuit.errors[aname] = self.errorParam[ipar]
            if self.paramStatusArray[ipar]==0:
                self.gMinuit.fixed[aname] = False
            else:
                self.gMinuit.fixed[aname] = True

        # set printlevel
        self.gMinuit.print_level = self.paramDict["printLevel"]

        # other setup
        self.nCallsCalcAll = 0


    def setupCoarseGrid(self,values):
        # set starting values of coarse grid
        self.startingParam = values
        self.errorParam = 0.01 * np.ones(self.npar)
        self.loParam = -1.0 * np.ones(self.npar) * self.paramDict["maxwavevalue"]
        self.hiParam = np.ones(self.npar) * self.paramDict["maxwavevalue"]
        self.paramStatusArray = np.zeros(self.npar)   # store =0 Floating, =1 Fixed

        # Set starting values and step sizes for parameters
        # (note that one can redefine the parameters, so this method can be called multiple times)
        self.parNames = []
        for ipar in range(self.npar):
            self.parNames.append("Grid_%d_%d" % (self.finegrid[ipar][0],self.finegrid[ipar][1]))


    def chisq(self,par):

        # call donutengine to calculate image

        # build coarse wavefront
        nCoarse = int(self._nbin/self._spacing)
        paramWave = np.zeros((nCoarse,nCoarse))
        for ipar in range(self.npar):
            j,i = self.coarsegrid[ipar]  #now consistent with assignment in this line: self.coarsegrid.append([j,i])
            paramWave[j,i] = par[ipar]

        # resize to full wavefront
        smoothWave = cv2.resize(paramWave,(self._nbin,self._nbin),interpolation=cv2.INTER_LANCZOS4)
        smoothWave = smoothWave * self.pupilMask

        # set in donutengine and calcImage!
        self.nCallsCalcAll +=1
        self.gFitFunc.setDeltaWFM(smoothWave)
        self.gFitFunc.calcAll(self.parCurrent)

        # compare to calculated image
        diff = self.imgarrayc - self.gFitFunc.getvImage()
        self.pullsq = diff*diff/self.sigmasq
        chisquared = self.pullsq.sum()

        # might need to add a counter-term to regularize the slopes of the paramWave grid?
        # or fit Fourier coefficients instead, to limit the the curve in frequency space...

        # printout
        if self.paramDict["printLevel"]>=2:
            # only print every 100 times
            if np.mod(self.nCallsCalcAll,100)==0:
                print('wavefit: Chi2 = ',chisquared," nCalls ",self.nCallsCalcAll)

        # save parameters for next iteration
        self._savePar = par.copy()
        self._chisq = chisquared

        # return result
        return chisquared

    def doFit(self):

        # start timer
        self.startingtime = time.time()

        # Migrad
        self.gMinuit.migrad(ncall=self.paramDict["maxIterations"])

        # done, check elapsed time
        firsttime = time.time()
        self.deltatime = firsttime - self.startingtime
        if self.paramDict["printLevel"]>=1:
            print('wavefit: Elapsed time fit = ',self.deltatime)

        # number of calls
        if self.paramDict["printLevel"]>=1:
            print('wavefit: Number of CalcAll calls = ',self.nCallsCalcAll)

    def outFit(self,postfix="wave",identifier=""):

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

        # print results
        print(tab.tabulate(*self.gMinuit.params.to_table()))

        # build coarse wavefront
        nCoarse = int(self._nbin/self._spacing)
        self.paramWave = np.zeros((nCoarse,nCoarse))
        self.paramWaveErr = np.zeros((nCoarse,nCoarse))
        for ipar in range(self.npar):
            j,i = self.coarsegrid[ipar]  # make this j,i; must be consistent!
            self.paramWave[j,i] = self.paramArray[ipar]
            self.paramWaveErr[j,i] = self.paramErrArray[ipar]


        # fill output Dictionary
        outputDict = {}
        outputDict["CHI2"] = amin
        outputDict["DOF"] = dof
        outputDict["FITSTAT"] = icstat
        outputDict["CLKTIME"] = self.deltatime
        outputDict["NCALLS"] = self.nCallsCalcAll

        for ipar in range(self.npar):
            outputDict[self.parNames[ipar]] = float(self.paramArray[ipar])
        for ipar in range(self.npar):
            outputDict[self.parNames[ipar]+"E"] = float(self.paramErrArray[ipar])

        # make a single output file, with multiple extensions
        # Extension 1:  Calculated Image
        # Extension 2:  Original Image
        # Extension 3:  Difference
        # Extension 4:  Chi2
        # Extension 5:  Wavefront

        hduListOutput = pyfits.HDUList()
        primaryOutput = pyfits.PrimaryHDU()
        primaryHeader = primaryOutput.header

        # fill primary header with fit results
        for key in outputDict:
            primaryHeader[key] = outputDict[key]
        hduListOutput.append(primaryOutput)

        # calculated Donut
        calcHdu = pyfits.ImageHDU(self.gFitFunc.getvImage(),name="Model_Image")
        hduListOutput.append(calcHdu)

        # original image
        imageHdu = pyfits.ImageHDU(self.imgarrayc,name="Original_Image")
        hduListOutput.append(imageHdu)

        # diff Donut - Calc
        diffHdu = pyfits.ImageHDU(self.imgarrayc-self.gFitFunc.getvImage(),name="Image-Model")
        hduListOutput.append(diffHdu)

        # Chi2 Donut-Calc
        chi2Hdu = pyfits.ImageHDU(self.pullsq,name="Chi2 Image-Model")
        hduListOutput.append(chi2Hdu)

        # Wavefront map - Zernike
        waveHdu = pyfits.ImageHDU(self.pupilMask*self.gFitFunc.getvPupilWaveZernike(),name="Zernike_Wavefront")
        hduListOutput.append(waveHdu)

        # Wavefront map - Plus Delta
        totalwavefront = self.pupilMask*self.gFitFunc.getvPupilWaveZernikePlusDelta()
        waveplusHdu = pyfits.ImageHDU(totalwavefront,name="Total_Wavefront")
        hduListOutput.append(waveplusHdu)

        # Wavefront map - Delta
        deltawavefront = self.pupilMask*(self.gFitFunc.getvPupilWaveZernikePlusDelta()-self.gFitFunc.getvPupilWaveZernike())
        wavedeltaHdu = pyfits.ImageHDU(deltawavefront,name="FineGrid_Wavefront")
        hduListOutput.append(wavedeltaHdu)

        # Wavefront map - Delta coarse
        wavedeltacoarseHdu = pyfits.ImageHDU(self.paramWave,name="CoarseGrid_Wavefront")
        hduListOutput.append(wavedeltacoarseHdu)

        # Wavefront map - Delta coarse errors
        wavedeltaerrcoarseHdu = pyfits.ImageHDU(self.paramWaveErr,name="Sigma_CoarseGrid_Wavefront")
        hduListOutput.append(wavedeltaerrcoarseHdu)

        # file names for output are
        # outputPrefix + identifier
        if identifier!="" :
            outName = self.paramDict["outputPrefix"] + "." + identifier + "." + postfix
        else:
            outName = self.paramDict["outputPrefix"] + "." + postfix
        outFile = outName + ".donut.fits"

        if self.paramDict["printLevel"]>=1:
            hduListOutput.info()
        hduListOutput.writeto(outFile,overwrite=True)

        # return some results
        resultDict = {}
        resultDict.update(outputDict)
        return resultDict
