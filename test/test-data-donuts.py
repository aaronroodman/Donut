###
### Fit one DECam donut with 3 method, first just Zernike polynomials, 
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

inprefix = 'input/DECam_00284692.S25.0020'
outprefix = 'output/DECam_00284692-fit%d.S25.0020'
infiles = [inprefix+'.stamp.fits']

# fit the donut three ways, just low order Zernikes, with wavefront map for high order Zernikes, and also with pupil grid

def fitdonut1():

    fitinitDict = {"nZernikeTerms":15,"fixedParamArray1":[0,1,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1],
                   "fixedParamArray2":[0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,0,0],"nFits":2,"nPixels":64,"nbin":512,
                   "scaleFactor":1.0,"pixelOverSample":8,"iTelescope":0,"inputrzero":0.15,"debugFlag":False,
                   "gain":4.5,"printLevel":1}

    df = donutfit(**fitinitDict)
    outfiles = [outprefix % (1)]

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
    outfiles = [outprefix % (2)]


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
    outfiles = [outprefix % (3)]

    # fit dictionary
    fitDict  = {}
    fitDict["inputrzero"] = 0.125
    fitDict["inputZernikeDict"] = {"S4":[0.0,0.0,11.0],"S25":[0.0,0.0,11.0],"None":[0.0,0.0,5.0]}

    for i in range(len(infiles)):
        fitDict["inputFile"] = infiles[i]
        fitDict["outputPrefix"] = outfiles[i]
        df.setupFit(**fitDict)

def anafits():
    outpostfix = '.second.donut.fits'
    files = [outprefix % (1)+outpostfix,outprefix % (2)+outpostfix,
outprefix % (3)+outpostfix]

    keys = ['CHI2','NELE','RZERO','BKGD']
    for iZ in range(4,15):
        keys.append('ZERN%d' % (iZ))

    outdict = {}
    for key in keys:
        outdict[key] = []

    for j,afile in enumerate(files):
        hdu = fits.open(afile)
        headd = hdu[0].header

        for i,key in enumerate(keys):
            outdict[key].append(headd[key])

    df = pd.DataFrame(outdict)
    print(tabulate(df, headers = 'keys',floatfmt='.3f'))


if __name__ == '__main__':
    fitdonut1()
    fitdonut2()
    fitdonut3()
    anafits()

