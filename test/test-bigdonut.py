###
### Script for fitting a BIG donut
###

import numpy as np
import scipy as sp
import pickle

from donutlib.donutfit import donutfit
from donutlib.wavefit import wavefit


fitinitDict = {"nZernikeTerms":15,"fixedParamArray1":[0,1,0,1,1,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1],"fixedParamArray2":[0,1,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1],"fixedParamArray3":[0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1],"nFits":3,"nPixels":256,"nbin":2048,"scaleFactor":1.0,"pixelOverSample":8,"iTelescope":0,"inputrzero":0.15,"outputWavefront":True,"debugFlag":False,"gain":4.5,"wavefrontMapFile" : "/Users/roodman/Astrophysics/Donuts/decam_2012-nominalzernike.pickle"}

df = donutfit(**fitinitDict)

# fit first donut
fitDict  = {}
fitDict["inputFile"] = 'DECam_00236392.S4.0003.stamp.fits'
#fitDict["inputFile"] = 'bigdonuts_00236392.S4.1.stamp.fits'
fitDict["outputPrefix"] = '/u/ec/roodman/kipacdisk/Donuts/bigdonuts/DECam_runwave_00236392.S4.0003'
fitDict["inputrzero"] = 0.125
fitDict["inputZernikeDict"] = {"S4":[0.0,0.0,53.0],"None":[0.0,0.0,11.0]}
df.setupFit(**fitDict)

df.gFitFunc.closeFits()

if True:
    # now fit an extra component of the wavefront, described by a mesh of points
    inputDict = {"outputPrefix":"wavetest","tolerance":3.0,"defineGrid":False,"spacing":64}

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

