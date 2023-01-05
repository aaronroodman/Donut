###
### Test making simulated Donuts with different wavelengths, but fitting with 700nm
###

from donutlib.makedonut import makedonut
from donutlib.donutfit import donutfit
import numpy as np
import pandas as pd
from tabulate import tabulate
from astropy.io import fits


def makefitw():

    # make donuts
    zarrnom = [0.,0.,0.,11.4,0.2,0.2,-0.15,0.05,0.3,-0.05,0.2,0.,0.,-0.07,0.10]

    # for DES
    inputDict = {'writeToFits':True,'outputPrefix':'unittest.wavelength.0001','iTelescope':0,'nZernikeTerms':15,'nbin':1024,'nPixels':64,'pixelOverSample':16,
                'scaleFactor':1.,'rzero':0.125, 'nEle':1.0e7, 'background':4000., 'randomFlag':False, 'randomSeed':2314809,
                'ZernikeArray':zarrnom,'gain':4.5,'printLevel':1,'waveLength':500.e-9}

    # build DataFrame
    headerTags = ['NELE','BKGD','RZERO','ZERN4','ZERN5','ZERN6','ZERN7','ZERN8','ZERN9','ZERN10','ZERN11','ZERN14','ZERN15',]
    dataf = pd.DataFrame(columns=headerTags)

    lams = [350,400,500,600,700,800,900,1000]
    for i,lam in enumerate(lams):
        inputDict["waveLength"] = lam * 1.e-9
        inputDict["outputPrefix"] = 'output/unittest.wavelength.%d' % (lam)
        inputDict['rzero'] = 0.125 * (lam/500.)   #is this correct?
        zarr = np.array(zarrnom) * (700./lam)
        inputDict['ZernikeArray'] = zarr.tolist()

        m = makedonut(**inputDict)
        donut1 = m.make()

        ## fit it

        fitinitDict = {"nZernikeTerms":15,"fixedParamArray1":[0,1,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1],
                    "fixedParamArray2":[0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,0,0],"nFits":2,"nPixels":64,"nbin":512,
                    "scaleFactor":1.0,"pixelOverSample":8,"iTelescope":0,"inputrzero":0.15,"debugFlag":False,"gain":4.5,"printLevel":1}

        df = donutfit(**fitinitDict)

        # fit donut
        fitDict  = {}
        fitDict["inputFile"] = 'output/unittest.wavelength.%d.stamp.fits' % (lam)
        fitDict["outputPrefix"] = 'output/unittest.wavelength.%d' % (lam)
        fitDict["inputrzero"] = 0.125
        fitDict["inputZernikeDict"] = {"N4":[0.0,0.0,5.0],"None":[0.0,0.0,11.0]}
        df.setupFit(**fitDict)

        # get fit Parameters
        hdu = fits.open("output/unittest.wavelength.%d.second.donut.fits" % (lam))
        header = hdu[0].header

        headList = []
        for aTag in headerTags:
            headList.append(header[aTag])

        dataf.loc[i] = headList

    print(tabulate(dataf, headers = 'keys',floatfmt='.3f'))

if __name__ == '__main__':
    makefitw()