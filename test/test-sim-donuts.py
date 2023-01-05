###
### Fit one simulated DECam donut with 3 methods, first just Zernike polynomials, 
### next with high order Zernikes from Zemax, last with pupil grid from BigDonuts
###

import os
import numpy as np
import pandas as pd
from astropy.io import fits
from tabulate import tabulate

from donutlib.makedonut import makedonut
from donutlib.donutfit import donutfit
from donutlib.wavefrontmap import wavefrontmap
from donutlib.decamutil import decaminfo

decamInfo = decaminfo()

wmap = os.path.expanduser("input/decam_2012-nominalzernike.pickle")
wmapobj = wavefrontmap(wmap)

def makeone():

    # make fake donut at ext='S4' and ix=500,iy=500
    xDECam,yDECam= decamInfo.getPosition(extname='S4',ix=500,iy=500)

    # get nominal wavefront 
    nZ = 37
    iZfirst = 4
    zarr = np.zeros(nZ+1)  # this array starts at iZ=0
    zarr[iZfirst:] = wmapobj.get(xDECam,yDECam,nZernikeFirst=iZfirst,nZernikeLast=nZ)

    # make a donut with slightly different Zernikes
    dzarr = np.zeros(nZ+1)
    dzarr[iZfirst:iZfirst+12] = 11.25,0.2,0.15,-0.15,0.05,0.3,-0.05,0.2,0.0,0.0,-0.025,0.025

    # make donut with pupil grid
    ztot = zarr + dzarr
    inputDict = {'xDECam':xDECam,'yDECam':yDECam,
        'writeToFits':True,'outputPrefix':'output/sim-donut.0001','iTelescope':0,'nZernikeTerms':37,'nbin':512,'nPixels':64,
        'pixelOverSample':8,'scaleFactor':1.,'rzero':0.20, 'nEle':1.0e7, 'background':500., 'randomFlag':True, 'randomSeed':2314809, 
        'ZernikeArray':ztot[1:],'gain':4.5,'printLevel':1,"deltaWFM":'input/DECam_236392_finegrid512.fits'}

    # make the donut
    m = makedonut(**inputDict)
    donut1 = m.make()

# fit the donut three ways, just low order Zernikes, with wavefront map for high order Zernikes, and also with pupil grid

def fitdonut1():

    fitinitDict = {"nZernikeTerms":15,"fixedParamArray1":[0,1,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1],
                   "fixedParamArray2":[0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,0,0],"nFits":2,"nPixels":64,"nbin":512,
                   "scaleFactor":1.0,"pixelOverSample":8,"iTelescope":0,"inputrzero":0.15,"debugFlag":False,
                   "gain":4.5,"printLevel":1}

    df = donutfit(**fitinitDict)
    infiles = ['output/sim-donut.0001.stamp.fits']
    outfiles = ['output/sim-donut-fit1.0001']

    # fit dictionary
    fitDict  = {}
    fitDict["inputrzero"] = 0.125
    fitDict["inputZernikeDict"] = {"S4":[0.0,0.0,11.0],"S25":[0.0,0.0,11.0],"None":[0.0,0.0,5.0]}

    for i in range(len(infiles)):
        fitDict["inputFile"] = infiles[i]
        fitDict["outputPrefix"] = outfiles[i]
        df.setupFit(**fitDict)

def fitdonut2():

    fitinitDict = {"nZernikeTerms":37,
        "fixedParamArray1":[0,1,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1],
        "fixedParamArray2":[0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1],
        "nFits":2,"nPixels":64,"nbin":512,
        "scaleFactor":1.0,"pixelOverSample":8,"iTelescope":0,"inputrzero":0.15,"debugFlag":False,
        "gain":4.5,"printLevel":1,"wavefrontMap":wmapobj}

    df = donutfit(**fitinitDict)
    infiles = ['output/sim-donut.0001.stamp.fits']
    outfiles = ['output/sim-donut-fit2.0001']

    # fit dictionary
    fitDict  = {}
    fitDict["inputrzero"] = 0.125
    fitDict["inputZernikeDict"] = {"S4":[0.0,0.0,11.0],"S25":[0.0,0.0,11.0],"None":[0.0,0.0,5.0]}

    for i in range(len(infiles)):
        fitDict["inputFile"] = infiles[i]
        fitDict["outputPrefix"] = outfiles[i]
        df.setupFit(**fitDict)

def fitdonut3():

    fitinitDict = {"nZernikeTerms":37,
        "fixedParamArray1":[0,1,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1],
        "fixedParamArray2":[0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1],
        "nFits":2,"nPixels":64,"nbin":512,
        "scaleFactor":1.0,"pixelOverSample":8,"iTelescope":0,"inputrzero":0.15,"debugFlag":False,
        "gain":4.5,"printLevel":1,"wavefrontMap":wmapobj,
        "deltaWFM":'input/DECam_236392_finegrid512.fits'}


    df = donutfit(**fitinitDict)
    infiles = ['output/sim-donut.0001.stamp.fits']
    outfiles = ['output/sim-donut-fit3.0001']

    # fit dictionary
    fitDict  = {}
    fitDict["inputrzero"] = 0.125
    fitDict["inputZernikeDict"] = {"S4":[0.0,0.0,11.0],"S25":[0.0,0.0,11.0],"None":[0.0,0.0,5.0]}

    for i in range(len(infiles)):
        fitDict["inputFile"] = infiles[i]
        fitDict["outputPrefix"] = outfiles[i]
        df.setupFit(**fitDict)

def anafits():
    files = ['output/sim-donut.0001.stamp.fits','output/sim-donut-fit1.0001.second.donut.fits','output/sim-donut-fit2.0001.second.donut.fits','output/sim-donut-fit3.0001.second.donut.fits']

    keysIN = ['NAXIS1','NELEINP','RZEROINP','BKGDINP']
    for iZ in range(4,15):
        keysIN.append('Z%d' % (iZ))

    keys = ['CHI2','NELE','RZERO','BKGD']
    for iZ in range(4,15):
        keys.append('ZERN%d' % (iZ))

    outdict = {}
    for key in keys:
        outdict[key] = []

    for j,afile in enumerate(files):
        hdu = fits.open(afile)
        headd = hdu[0].header
        if j==0:
            for i,key in enumerate(keys):
                keyIN = keysIN[i]
                outdict[key].append(headd[keyIN])
        else:
            for i,key in enumerate(keys):
                outdict[key].append(headd[key])

    df = pd.DataFrame(outdict)
    print(tabulate(df, headers = 'keys',floatfmt='.3f'))


if __name__ == '__main__':
    makeone()
    fitdonut1()
    fitdonut2()
    fitdonut3()
    anafits()

