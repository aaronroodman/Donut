###
### Script for very basic testing of donutengine and donutfit
###

from donutlib.makedonut import makedonut
from donutlib.donutfit import donutfit
from donutlib.wavefrontmap import wavefrontmap
import pdb
# make a donut

# zernike coeff at sensor N1 (change z4 by +11.0)
zarray = [0.16335417, 0.04209118, 0.47722926, 11.0+0.01731888, -0.02072751, -0.11661182, -0.07499869, 0.00662177, 0.26644876, -0.07200162, -0.0686725, 0.04443583, 0.00789245, 0.03902099, 0.01432855, -0.00331504, 0.0375199, -0.01361314, 0.0503574, 0.00222714, -0.00473283, -0.00852732, -0.00527068, -0.02972607, 0.00167757, 0.0045495, -0.00015644, -0.00026905, -0.00856108, 0.00075652, 0.00293147, -0.00078142, -0.00017493, 8.541e-05, 6.06e-06, -3.87e-06, 0.0013547]

xdecam,ydecam = 16.9203, -191.826

wmap = "/u/ec/roodman/Astrophysics/Donuts/decam_2012-nominalzernike.pickle"
wmapobj = wavefrontmap(wmap)

inputDict = {'writeToFits':True,'outputPrefix':'unittestw.0001','iTelescope':0,'nZernikeTerms':37,'nbin':512,'nPixels':64,'pixelOverSample':8,'scaleFactor':1.,'rzero':0.125, 'nEle':1.0e6, 'background':4000., 'randomFlag':True, 'randomSeed':2314809, 'ZernikeArray':zarray,'gain':4.5,'xDECam':xdecam,'yDECam':ydecam}

m = makedonut(**inputDict)
donut1 = m.make()

## fit them

fitinitDict = {"nZernikeTerms":37,"fixedParamArray1":[0,1,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1],"fixedParamArray2":[0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1],"nFits":2,"nPixels":64,"nbin":512,"scaleFactor":1.0,"pixelOverSample":8,"iTelescope":0,"inputrzero":0.15,"debugFlag":False,"gain":4.5,"wavefrontMap":wmapobj}

df = donutfit(**fitinitDict)

# fit first donut
fitDict  = {}
fitDict["inputFile"] = 'unittestw.0001.stamp.fits'
fitDict["outputPrefix"] = 'unittestw.0001'
fitDict["inputrzero"] = 0.125
fitDict["inputZernikeDict"] = {"N1":[0.0,0.0,11.0],"None":[0.0,0.0,11.0]}
df.setupFit(**fitDict)


