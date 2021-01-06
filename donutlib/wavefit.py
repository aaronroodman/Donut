import time
import os
import numpy
import sys
import cv2
from collections import OrderedDict
from astropy.io import fits as pyfits
from array import array
import pdb
import ctypes

import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True
from donutlib.donutengine import donutengine

class wavefit(object):
    """ donutfit is a class used to fit the wavefront of donuts, using donutengine and MINUIT, for the DES experiment

    Aaron Roodman (C) SLAC National Accelerator Laboratory, Stanford University 2019.
    """

    def __init__(self,dfObject,**inputDict):
        # init contains all initializations which are done only once for all fits
        # parameters in fixParamArray1 are nEle,rzero,bkgd,Z2,Z3,Z4,....Z11
        
        self.paramDict = {"outputPrefix":"test","printLevel":2,"maxIterations":100,"spacing":16,"chi2factor1":1.e5,"maxwavevalue":0.5,"tolerance":0.1,"defineGrid":True}

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
                 
        # get weights and background subtracted image 
        constantError2 = 7.1 * 7.1 
        self.sigmasq = self.df.imgarray + constantError2
        self.weight = 1.0/self.sigmasq

        self.imgarrayc = self.df.imgarray - bkg        

        # setup MINUIT
        self.gMinuit = ROOT.TMinuit(self.npar)
        self.gMinuit.SetFCN( self.chisq )

        # arglist is for the parameters in Minuit commands
        arglist = array( 'd', 10*[0.] )
        ierflg =ctypes.c_int(1982)#L ROOT.Long(1982) 

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

        if self.paramDict["defineGrid"]:
            # status/limit arrays for Minuit parameters
            self.startingParam = numpy.zeros(self.npar)
            self.errorParam = 0.01 * numpy.ones(self.npar)
            self.loParam = -1.0 * numpy.ones(self.npar) * self.paramDict["maxwavevalue"]
            self.hiParam = numpy.ones(self.npar) * self.paramDict["maxwavevalue"]
            self.paramStatusArray = numpy.zeros(self.npar)   # store =0 Floating, =1 Fixed

            # Set starting values and step sizes for parameters
            # (note that one can redefine the parameters, so this method can be called multiple times)
            self.parNames = []
            for ipar in range(self.npar):
                self.parNames.append("Grid_%d_%d" % (self.finegrid[ipar][0],self.finegrid[ipar][1]))
                self.gMinuit.DefineParameter(ipar,"Grid_%d_%d" % (self.finegrid[ipar][0],self.finegrid[ipar][1]),self.startingParam[ipar],self.errorParam[ipar],self.loParam[ipar],self.hiParam[ipar])

        # other setup
        self.nCallsCalcAll = 0


    def setupCoarseGrid(self,values):
        # set starting values of coarse grid
        self.startingParam = values
        self.errorParam = 0.01 * numpy.ones(self.npar)
        self.loParam = -1.0 * numpy.ones(self.npar) * self.paramDict["maxwavevalue"]
        self.hiParam = numpy.ones(self.npar) * self.paramDict["maxwavevalue"]
        self.paramStatusArray = numpy.zeros(self.npar)   # store =0 Floating, =1 Fixed

        # Set starting values and step sizes for parameters
        # (note that one can redefine the parameters, so this method can be called multiple times)
        self.parNames = []
        for ipar in range(self.npar):
            self.parNames.append("Grid_%d_%d" % (self.finegrid[ipar][0],self.finegrid[ipar][1]))
            self.gMinuit.DefineParameter(ipar,"Grid_%d_%d" % (self.finegrid[ipar][0],self.finegrid[ipar][1]),self.startingParam[ipar],self.errorParam[ipar],self.loParam[ipar],self.hiParam[ipar])



                
    def chisq(self,npar, gin, f, par, iflag ):

        # convert par to a numpy array
        parArr = numpy.zeros(self.npar)
        for ipar in range(self.npar):
            parArr[ipar] = par[ipar]
            
        # call donutengine to calculate image

        # build coarse wavefront
        nCoarse = int(self._nbin/self._spacing)
        paramWave = numpy.zeros((nCoarse,nCoarse))
        for ipar in range(self.npar):
            i,j = self.coarsegrid[ipar]
            paramWave[j,i] = parArr[ipar] 

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
            print('donutfit: Chi2 = ',chisquared)

        # save parameters for next iteration
        self._savePar = parArr
        self._chisq = chisquared

        # return result    
        f[0] = chisquared

    def doFit(self):
        
        # arglist is for the parameters in Minuit commands
        arglist = array( 'd', 10*[0.] )
        ierflg = ctypes.c_int(1982) #L ROOT.Long

        # tell Minuit to use strategy for best fits
        arglist[0] = 1  # 
        self.gMinuit.mnexcm( "SET STRATEGY", arglist, 1, ierflg )
                
        # start timer
        self.startingtime = time.clock()

        # Now ready for minimization step
        #self.gMinuit.SetMaxIterations(self.paramDict["maxIterations"])
        #self.gMinuit.Migrad()

        # Migrad
        arglist[0] = 10000  # maxcalls
        arglist[1] = self.paramDict["tolerance"]    # tolerance, default is 0.1
        self.gMinuit.mnexcm( "MIGRAD", arglist, 2, ierflg )
        

        # done, check elapsed time
        firsttime = time.clock()
        self.deltatime = firsttime - self.startingtime
        if self.paramDict["printLevel"]>=1:
            print('wavefit: Elapsed time fit = ',self.deltatime)

        # number of calls
        if self.paramDict["printLevel"]>=1:
            print('wavefit: Number of CalcAll calls = ',self.nCallsCalcAll)

    def outFit(self,postfix=".wave.donut.fits",identifier=""):

        # get more fit details from MINUIT
        amin, edm, errdef = ctypes.c_double(0.18), ctypes.c_double(0.19), ctypes.c_double(0.20)
        nvpar, nparx, icstat = ctypes.c_int(1983), ctypes.c_int(1984), ctypes.c_int(1985)
        self.gMinuit.mnstat( amin, edm, errdef, nvpar, nparx, icstat )
        dof = pow(self._nPixels,2) - nvpar.value
        if self.paramDict["printLevel"]>=1:
            mytxt = "amin = %.3f, edm = %.3f,   effdef = %.3f,   nvpar = %.3f,  nparx = %.3f, icstat = %.3f " % (amin.value,edm.value,errdef.value,nvpar.value,nparx.value,icstat.value)   
            print('wavefit: ',mytxt)

        # get fit values and errors
        aVal = ctypes.c_double(0.21)
        errVal = ctypes.c_double(0.22)
        self.paramArray = numpy.zeros(self.npar)
        self.paramErrArray = numpy.zeros(self.npar)
        for ipar in range(self.npar):
            self.gMinuit.GetParameter(ipar,aVal,errVal)
            self.paramArray[ipar] = aVal.value
            if errVal.value < 1e9 :
                self.paramErrArray[ipar] = errVal.value
            else:
                self.paramErrArray[ipar] = 0.0

        # build coarse wavefront
        nCoarse = int(self._nbin/self._spacing)
        self.paramWave = numpy.zeros((nCoarse,nCoarse))
        self.paramWaveErr = numpy.zeros((nCoarse,nCoarse))
        for ipar in range(self.npar):
            i,j = self.coarsegrid[ipar]
            self.paramWave[j,i] = self.paramArray[ipar]
            self.paramWaveErr[j,i] = self.paramErrArray[ipar]
                        
                
        # fill output Dictionary
        outputDict = {}
        outputDict["CHI2"] = float(amin.value)
        outputDict["DOF"] = int(dof)
        outputDict["FITSTAT"] = int(icstat.value)
        outputDict["CLKTIME"] = self.deltatime
        outputDict["DOF"] = dof
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
        
        # file names for output 
        outName = self.paramDict["outputPrefix"] 

        # write out fits file
        outFile =  outName + postfix
        if self.paramDict["printLevel"]>=1:
            hduListOutput.info()
        hduListOutput.writeto(outFile,overwrite=True)
        
        # return some results
        resultDict = {}
        resultDict.update(outputDict)
        return resultDict




