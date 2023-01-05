###
### Fit 4 donuts with just Zernike polynomials
###

from donutlib.donutfit import donutfit

def fitdonut():

    fitinitDict = {"nZernikeTerms":15,"fixedParamArray1":[0,1,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1],
                   "fixedParamArray2":[0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,0,0],"nFits":2,"nPixels":64,"nbin":512,
                   "scaleFactor":1.0,"pixelOverSample":8,"iTelescope":0,"inputrzero":0.15,"debugFlag":False,
                   "gain":4.5,"printLevel":1}

    df = donutfit(**fitinitDict)

    infiles = ['input/DECam_00284696.S4.0008.stamp.fits','input/DECam_00284696.S25.0001.stamp.fits',
               'input/DECam_00345461.S4.0018.stamp.fits','input/DECam_00345461.S25.0018.stamp.fits']

    outfiles = ['output/DECam_00284696.S4.0008','output/DECam_00284696.S25.0001',
                'output/DECam_00345461.S4.0018','output/DECam_00345461.S25.0018']

    # fit dictionary
    fitDict  = {}
    fitDict["inputrzero"] = 0.125
    fitDict["inputZernikeDict"] = {"S4":[0.0,0.0,11.0],"S25":[0.0,0.0,11.0],"None":[0.0,0.0,5.0]}

    for i in range(len(infiles)):
        fitDict["inputFile"] = infiles[i]
        fitDict["outputPrefix"] = outfiles[i]
        df.setupFit(**fitDict)


if __name__ == '__main__':
    fitdonut()