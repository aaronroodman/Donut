#
# test pupil grid fit on a regular donut, save pupil grid and also debug DonutEngine output
#
import os
import numpy as np
import scipy as sp
import pickle

from donutlib.donutfit import donutfit
from donutlib.wavefit import wavefit
from donutlib.wavefrontmap import wavefrontmap

wmap = os.path.expanduser("input/decam_2012-nominalzernike.pickle")
wmapobj = wavefrontmap(wmap)

def fitsmalldonut():

    fitinitDict = {"nZernikeTerms":37,"fixedParamArray1":[0,1,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1],
                "fixedParamArray2":[0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1],
                "nFits":2,"nPixels":64,"nbin":512,"scaleFactor":1.0,"pixelOverSample":8,"iTelescope":0,"inputrzero":0.15,"outputWavefront":True,"debugFlag":True,"gain":4.5,"wavefrontMap":wmapobj}

    df = donutfit(**fitinitDict)

    # fit first donut
    fitDict  = {}
    fitDict["inputFile"] = 'input/DECam_00284696.S4.0008.stamp.fits'
    fitDict["outputPrefix"] = 'output/DECam_pupilgridfit_00284696.S4.0008'
    fitDict["inputrzero"] = 0.125
    fitDict["inputZernikeDict"] = {"S4":[0.0,0.0,11.0],"None":[0.0,0.0,11.0]}
    df.setupFit(**fitDict)

    df.gFitFunc.closeFits()

    # now fit an extra component of the wavefront, described by a mesh of points
    inputDict = {"outputPrefix":fitDict["outputPrefix"],"maxIterations":30000,"tolerance":400000.0}

    wfit = wavefit(df,**inputDict)

    # do the fit
    wfit.doFit()

    # get the results
    wfit.outFit()


if __name__ == '__main__':
    fitsmalldonut()   