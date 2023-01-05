###
### Make and Fit one simulated donut
###

from donutlib.makedonut import makedonut
from donutlib.donutfit import donutfit

def makeandfit():

    # make donuts
    z4 =11.
    z5 = 0.2
    z6 = 0.2
    z7 = -0.15
    z8 = .05
    z9 = 0.3
    z10 = -.05
    z11 = .2


    inputDict = {'writeToFits':True,'outputPrefix':'output/test-sim.0001','iTelescope':0,
                 'nZernikeTerms':37,'nbin':512,'nPixels':64,'pixelOverSample':8,'scaleFactor':1.,'rzero':0.125, 
                 'nEle':1.0e6, 'background':4000., 'randomFlag':True, 'randomSeed':2314809, 
                 'ZernikeArray':[0.,0.,0.,z4,z5,z6,z7,z8,z9,z10,z11],'gain':4.5,'printLevel':1}


    m = makedonut(**inputDict)
    donut1 = m.make()


    ## fit it

    fitinitDict = {"nZernikeTerms":15,"fixedParamArray1":[0,1,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1],
                   "fixedParamArray2":[0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,0,0],"nFits":2,"nPixels":64,"nbin":512,
                   "scaleFactor":1.0,"pixelOverSample":8,"iTelescope":0,"inputrzero":0.15,"debugFlag":True,
                   "gain":4.5,"printLevel":1}

    df = donutfit(**fitinitDict)

    # fit first donut
    fitDict  = {}
    fitDict["inputFile"] = 'output/test-sim.0001.stamp.fits'
    fitDict["outputPrefix"] = 'output/test-sim.0001'
    fitDict["inputrzero"] = 0.125
    fitDict["inputZernikeDict"] = {"N4":[0.0,0.0,11.0],"None":[0.0,0.0,11.0]}
    df.setupFit(**fitDict)


if __name__ == '__main__':
    makeandfit()