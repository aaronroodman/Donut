###
### Compare fits of 4 Donuts with 3 methods: low order Zernikes, with Zemax high order Zernikes, and also with pupil grid
###

import os
import numpy as np
import pandas as pd
from astropy.io import fits
from tabulate import tabulate

filenames = ['output/DECam%s_00284696.S25.0001.second.donut.fits',
                'output/DECam%s_00284696.S4.0008.second.donut.fits',
                'output/DECam%s_00345461.S25.0018.second.donut.fits',
                'output/DECam%s_00345461.S4.0018.second.donut.fits']

fitnames = ["","_wmap","_pupilgrid"]

def anafits():

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
    anafits()

