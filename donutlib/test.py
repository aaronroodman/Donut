###
### Script for very basic testing of donutengine and donutfit
###

from donutlib.makedonut import makedonut
from donutlib.donutfit import donutfit


# make donuts
z4 = 10.
z5 = 0.2
z6 = 0.2
z7 = -0.15
z8 = .05
z9 = 0.3
z10 = -.05
z11 = .2
inputDict = {'writeToFits':True,'outputPrefix':'unittest.0001','iTelescope':5,'nZernikeTerms':37,'nbin':512,'nPixels':64,'pixelOverSample':8,'scaleFactor':1.,'rzero':0.125, 'nEle':1.0e6, 'background':4000., 'randomFlag':True, 'randomSeed':2314809, 'ZernikeArray':[0.,0.,0.,z4,z5,z6,z7,z8,z9,z10,z11],'xDECam':1.318,'yDECam':0.86}

m = makedonut(**inputDict)
donut1 = m.make()

z4 = 10.5
z9 = -0.05
z10 = 0.3
newDict = {'outputPrefix':'unittest.0002','ZernikeArray':[0.,0.,0.,z4,z5,z6,z7,z8,z9,z10,z11]}
donut2 = m.make(**newDict)

## fit them
fitinitDict = {"nZernikeTerms":15,"fixedParamArray1":[0,1,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1],"fixedParamArray2":[0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,0,0],"nFits":2,"nPixels":64,"nbin":512,"scaleFactor":1.0,"pixelOverSample":8,"iTelescope":5,"inputrzero":0.15,"debugFlag":False}
df = donutfit(**fitinitDict)

# fit first donut
fitDict  = {}
fitDict["inputFile"] = 'unittest.0001.stamp.fits'
fitDict["outputPrefix"] = 'unittest.0001'
fitDict["inputrzero"] = 0.125
fitDict["inputZernikeDict"] = {"N4":[0.0,0.0,11.0]}
df.setupFit(**fitDict)

# fit second one
fitDict["inputFile"] = 'unittest.0002.stamp.fits'
fitDict["outputPrefix"] = 'unittest.0002'
df.setupFit(**fitDict)



