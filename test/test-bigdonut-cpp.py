###
### Script for fitting a BIG donut
###

import numpy as np

from donutlib.donutfit import donutfit

fitinitDict = {"nZernikeTerms":15,"fixedParamArray1":[0,1,0,1,1,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1],"fixedParamArray2":[0,1,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1],"fixedParamArray3":[0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1],"nFits":3,"nPixels":256,"nbin":2048,"scaleFactor":1.0,"pixelOverSample":8,"iTelescope":0,"inputrzero":0.15,"outputWavefront":True,"debugFlag":False,"gain":4.5,"wavefrontMapFile" : "/Users/roodman/Astrophysics/Donuts/decam_2012-nominalzernike.pickle", "doGridFit":True, "spacing":64}

df = donutfit(**fitinitDict)

# fit donut
fitDict  = {}
fitDict["inputFile"] = 'DECam_00236392.S4.0003.stamp.fits'
fitDict["outputPrefix"] = 'DECam_wave_00236392.S4.0003'
fitDict["inputrzero"] = 0.125
fitDict["inputZernikeDict"] = {"S4":[0.0,0.0,53.0],"None":[0.0,0.0,11.0]}
df.setupFit(**fitDict)

df.gFitFunc.closeFits()

