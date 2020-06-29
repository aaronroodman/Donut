
import matplotlib
matplotlib.rcParams['agg.path.chunksize'] = 10000
from matplotlib import pyplot as plt

import numpy as np
import scipy as sp
import pickle

from donutlib.donutfit import donutfit
from donutlib.wavefit import wavefit


fitinitDict = {"nZernikeTerms":15,"fixedParamArray1":[0,1,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1],"fixedParamArray2":[0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1],"nFits":2,"nPixels":64,"nbin":512,"scaleFactor":1.0,"pixelOverSample":8,"iTelescope":0,"inputrzero":0.15,"outputWavefront":True,"debugFlag":True,"gain":4.5,"wavefrontMapFile" : "/Users/roodman/Astrophysics/Donuts/decam_2012-nominalzernike.pickle"}

df = donutfit(**fitinitDict)

# fit first donut
fitDict  = {}
fitDict["inputFile"] = 'DECam_00284696.S4.0008.stamp.fits'
fitDict["outputPrefix"] = 'DECam_runwave_00284696.S4.0008'
fitDict["inputrzero"] = 0.125
fitDict["inputZernikeDict"] = {"S4":[0.0,0.0,11.0],"None":[0.0,0.0,11.0]}
df.setupFit(**fitDict)

df.gFitFunc.closeFits()

# now fit an extra component of the wavefront, described by a mesh of points
inputDict = {"outputPrefix":"wavetest","tolerance":3.0,"defineGrid":False}

wfit = wavefit(df,**inputDict)

# setup initial values
values = np.zeros(wfit.npar)
for ipar in range(wfit.npar):
    ix,iy = wfit.coarsegrid[ipar]  # in case we want starting values to vary with x,y
    values[ipar] = 0.0
wfit.setupCoarseGrid(values)

# do the fit
wfit.doFit()

# get the results
wfit.outFit()


#----

#        hduListOutput = pyfits.HDUList()
#        primaryOutput = pyfits.PrimaryHDU()
#        primaryHeader = primaryOutput.header

        # fill primary header with fit results
#        for key in outputDict:
#            primaryHeader[key] = outputDict[key]
#        hduListOutput.append(primaryOutput)
        
        # calculated Donut
#        calcHdu = pyfits.ImageHDU(wfit.gFitFunc.getvImage())
#        hduListOutput.append(calcHdu)

        # original image
#        imageHdu = pyfits.ImageHDU(wfit.imgarrayc)
#        hduListOutput.append(imageHdu)
            
        # diff Donut - Calc
#        diffHdu = pyfits.ImageHDU(wfit.imgarrayc-wfit.gFitFunc.getvImage())
#        hduListOutput.append(diffHdu)

        # Chi2 Donut-Calc
#        chi2Hdu = pyfits.ImageHDU(wfit.pullsq)
#        hduListOutput.append(chi2Hdu)

        # Wavefront map - Zernike
#        waveHdu = pyfits.ImageHDU(wfit.pupilMask*wfit.gFitFunc.getvPupilWaveZernike())
#        hduListOutput.append(waveHdu)

        # Wavefront map - Plus Delta
#        waveplusHdu = pyfits.ImageHDU(wfit.pupilMask*wfit.gFitFunc.getvPupilWaveZernikePlusDelta())
#        hduListOutput.append(waveplusHdu)

        # Wavefront map - Delta
#        wavedeltaHdu = pyfits.ImageHDU(wfit.pupilMask*(wfit.gFitFunc.getvPupilWaveZernikePlusDelta()-wfit.gFitFunc.getvPupilWaveZernike()))
#        hduListOutput.append(wavedeltaHdu)


# make star with and without

# plt.interactive(True)

# gfit = wfit.gFitFunc
# par = gfit.getParCurrent()
# par[5] = .1
# gfit.calcAll(par)
# star_withdelta = gfit.getvImage()

# f = plt.figure()
# plt.imshow(star_withdelta,origin='lower',interpolation='None')

# gfit.unsetDeltaWFM()

# zero_dwfm = np.zeros((512,512))
# gfit.setDeltaWFM(zero_dwfm)

# gfit.calcAll(par)
# star_withoutdelta = gfit.getvImage()

# f = plt.figure()
# plt.imshow(star_withoutdelta,origin='lower',interpolation='None')

# diff = (star_withdelta - star_withoutdelta)/np.sum(star_withoutdelta)
# f = plt.figure()
# plt.imshow(diff,origin='lower',interpolation='None')


# pwzpd = gfit.getvPupilWaveZernikePlusDelta()
# pwz = gfit.getvPupilWaveZernike()
# pd = pwz - pwzpd

