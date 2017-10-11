""" 24 August 2017. Author: Sasha Safonova.
Script that directly calls the calcPupilMask method in DonutEngine.cc, created for development purposes only.
Be sure to compile donutengine first by:

cd ../src
make clean
make swig
make
"""


from donutlib.donutengine import donutengine
import numpy as np
from matplotlib import pyplot as plt
from skimage.transform import resize
import pandas as pd


def readimg(filename, nbin):
    """Reads an example wavefront stamp in the "trim" format and turns it into a flat pupil image.

    Parameters
    ----------
    filename : string
        full path to the text file containing the wavefront stamp 
    nbin : int
        the number of pixels that the output stamp should have per side 


    Returns
    -------
    img : array
        the final pupil stamp with the number of pixels per side defined by nbin
    """
    img = pd.read_csv(filename, sep='\s+', header=None)
    img = img.as_matrix()
    img[img != 0] = 1
    img = resize(img, (nbin, nbin), mode='reflect')
    return img

# Set the number of pixels per side in the stamp, e.g. 512, 256, 128, etc.
binns = 512

# Set up a donutengine object
paramDict =  {"inputFile":"",
                          "wfmFile":"",
                          "wfmArray":None,
                          "writeToFits":False,
                          "outputPrefix":"testone",
                          "xDECam":0,
                          "yDECam":0,
                          "debugFlag":False,
                          "rootFlag":False,
                          "waveLength":700.0e-9,
                          "nZernikeTerms":37,
                          "nbin":binns,
                          "nPixels":64,
                          "gridCalcMode":True,
                          "pixelOverSample":8,
                          "scaleFactor":1.,                 
                          "rzero":0.125,
                          "nEle":1.0e6,
                          "background":4000.,
                          "randomFlag":False,
                          "randomSeed":209823,  # if this is an invalid integer, crazy errors will ensue
                          "gain":1.0,
                          "flipFlag":False,
                          "iTelescope": 5,      # 5 stands for the DESI configuration
                          "ZernikeArray":[]}
gFitFunc = donutengine(**paramDict)
parArr = np.zeros(gFitFunc.npar)
parArr[gFitFunc.ipar_bkgd] = 4000.
parArr[gFitFunc.ipar_nEle] = 1.e6
parArr[gFitFunc.ipar_rzero] = 0.15
parArr[gFitFunc.ipar_ZernikeFirst+2] = 7.4

# Set the angle (in degrees, Zemax Echo 22 coordinate system)
# See angles for the test fields inside 00README in the testpupils directory
gFitFunc.setXYDESI(1.325, 0.89)

#calculate and get the pupil mask
gFitFunc.calcAll(parArr)
data = gFitFunc.getvImage()
pupil = gFitFunc.getvPupilMask()

# Change the test wavefront file to be used for comparison here
filename = './testpupils/field11_trim.txt'
comparisondonut=readimg(filename, binns)

plt.figure()
#plot the donutengine output in color and the ray-tracing donut over it in grey
plt.imshow(pupil, alpha=0.5)  
plt.imshow(comparisondonut, alpha=0.5, cmap='Greys')
plt.ylim(0, binns)
plt.savefig("seetheplot.png")
plt.show()

