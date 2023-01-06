###
### Fit 4 donuts with Zernike polynomials and high order Zernikes from Zemax and pupil grid from Big Donuts
###

import os
import numpy as np
import pandas as pd
from astropy.io import fits
from tabulate import tabulate
from donutlib.donutfit import donutfit
from donutlib.wavefrontmap import wavefrontmap


infiles = ['input/DECam_00284696.S4.0008.stamp.fits','input/DECam_00284696.S25.0001.stamp.fits',
           'input/DECam_00345461.S4.0018.stamp.fits','input/DECam_00345461.S25.0018.stamp.fits',
            'input/DECam_00153648.S4.0010.stamp.fits','input/DECam_00153648.S25.0019.stamp.fits']

filenames = ['out'+afile[2:11]+'%s'+afile[11:-11]+'.second.donut.fits' for afile in infiles]

def fitdonut1():

    fitinitDict = {"nZernikeTerms":15,"fixedParamArray1":[0,1,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1],
                   "fixedParamArray2":[0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,0,0],"nFits":2,"nPixels":64,"nbin":512,
                   "scaleFactor":1.0,"pixelOverSample":8,"iTelescope":0,"inputrzero":0.15,"debugFlag":False,
                   "gain":4.5,"printLevel":1}

    df = donutfit(**fitinitDict)

    outfiles = ['out'+afile[2:-11] for afile in infiles]

    # fit dictionary
    fitDict  = {}
    fitDict["inputrzero"] = 0.125
    fitDict["inputZernikeDict"] = {"S4":[0.0,0.0,11.0],"S25":[0.0,0.0,11.0],"None":[0.0,0.0,5.0]}

    for i in range(len(infiles)):
        fitDict["inputFile"] = infiles[i]
        fitDict["outputPrefix"] = outfiles[i]
        df.setupFit(**fitDict)

def fitdonut2():

    wmap = os.path.expanduser("input/decam_2012-nominalzernike.pickle")
    wmapobj = wavefrontmap(wmap)

    fitinitDict = {"nZernikeTerms":37,
        "fixedParamArray1":[0,1,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1],
        "fixedParamArray2":[0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1],
        "nFits":2,"nPixels":64,"nbin":512,
        "scaleFactor":1.0,"pixelOverSample":8,"iTelescope":0,"inputrzero":0.15,"debugFlag":False,
        "gain":4.5,"printLevel":1,"wavefrontMap":wmapobj}

    df = donutfit(**fitinitDict)

    outfiles = ['out'+afile[2:11]+'_wmap'+afile[11:-11] for afile in infiles]

    # fit dictionary
    fitDict  = {}
    fitDict["inputrzero"] = 0.125
    fitDict["inputZernikeDict"] = {"S4":[0.0,0.0,11.0],"S25":[0.0,0.0,11.0],"None":[0.0,0.0,5.0]}

    for i in range(len(infiles)):
        fitDict["inputFile"] = infiles[i]
        fitDict["outputPrefix"] = outfiles[i]
        df.setupFit(**fitDict)        

def fitdonut3():

    wmap = os.path.expanduser("input/decam_2012-nominalzernike.pickle")
    wmapobj = wavefrontmap(wmap)

    fitinitDict = {"nZernikeTerms":37,
        "fixedParamArray1":[0,1,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1],
        "fixedParamArray2":[0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1],
        "nFits":2,"nPixels":64,"nbin":512,
        "scaleFactor":1.0,"pixelOverSample":8,"iTelescope":0,"inputrzero":0.15,"debugFlag":False,
        "gain":4.5,"printLevel":1,"wavefrontMap":wmapobj,
        "deltaWFM":'input/DECam_236392_finegrid512.fits'}

    df = donutfit(**fitinitDict)

    outfiles = ['out'+afile[2:11]+'_pupilgrid'+afile[11:-11] for afile in infiles]

    # fit dictionary
    fitDict  = {}
    fitDict["inputrzero"] = 0.125
    fitDict["inputZernikeDict"] = {"S4":[0.0,0.0,11.0],"S25":[0.0,0.0,11.0],"None":[0.0,0.0,5.0]}

    for i in range(len(infiles)):
        fitDict["inputFile"] = infiles[i]
        fitDict["outputPrefix"] = outfiles[i]
        df.setupFit(**fitDict)


def anafits():

    fitnames = ["","_wmap","_pupilgrid"]

    keys = ['CHI2','NELE','RZERO','BKGD']
    for iZ in range(4,15):
        keys.append('ZERN%d' % (iZ))

    outdict = {}
    outdict['FITNAME'] = []
    outdict['FILENAME'] = []

    for key in keys:
        outdict[key] = []

    for j,afile in enumerate(filenames):
        for k,afit in enumerate(fitnames):

            thefile = afile % (afit)
            hdu = fits.open(thefile)
            headd = hdu[0].header
            outdict['FITNAME'].append(afit)
            outdict['FILENAME'].append(afile[15:33])

            for i,key in enumerate(keys):
                outdict[key].append(headd[key])

    df = pd.DataFrame(outdict)
    print(tabulate(df, headers = 'keys',floatfmt='.3f'))


if __name__ == '__main__':
    fitdonut1()
    fitdonut2()
    fitdonut3()
    anafits()