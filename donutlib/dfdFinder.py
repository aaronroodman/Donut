import numpy
import time
from numpy.lib.stride_tricks import as_strided as ast
from scipy.spatial import cKDTree
from donutlib.decamutil import decaminfo
import pdb

def block_view(A, block= (8, 8)):
    """Provide a 2D block view to 2D array. No error checking made.
    Therefore meaningful (as implemented) only for blocks strictly
    compatible with the shape of A."""
    # simple shape and strides computations may seem at first strange
    # unless one is able to recognize the 'tuple additions' involved ;-)
    shape= (A.shape[0]// block[0], A.shape[1]// block[1])+ block
    strides= (block[0]* A.strides[0], block[1]* A.strides[1])+ A.strides
    return ast(A, shape= shape, strides= strides)


def calcSky(imgarray,nsigma=1.5,debugFlag=False):

    # size 
    nArea = imgarray.shape[0]*imgarray.shape[1]

    # get Mean,Sigma, Mask of pixels below Mean+nsigma*sigma
    # use numpy Masked Array!
    imgMean = imgarray.mean()
    imgStd = imgarray.std()
    if debugFlag:
        print("calcSky: mean,sigma = ",imgMean,imgStd)
    imgarrayMask = numpy.ma.masked_greater(imgarray,imgMean+1*imgStd)
    countsNew = (imgarrayMask.mask==False).sum()
    countsOld = -1

    # only try 4 times at maximum!
    itry = 0

    while countsOld!=countsNew and itry<4 :
        itry = itry + 1
        countsOld = countsNew
        imgMean = imgarrayMask.mean()
        imgStd = imgarrayMask.std()
        if debugFlag:
            print("calcSky: mean,sigma,counts = ",imgMean,imgStd,countsOld)
        imgarrayMask = numpy.ma.masked_greater(imgarray,imgMean+nsigma*imgStd)
        countsNew = (imgarrayMask.mask==False).sum()
                
    return imgMean


def findDonuts(dataAmp,xOffset,inputDict,extName):
    """ find donuts in the data array
    """

    nBlock = inputDict["nBlock"]
    fluxThreshold = inputDict["fluxThreshold"]
    kNNFind = inputDict["kNNFind"]
    distanceLimit = inputDict["distanceLimit"]
    fluxMinCut = inputDict["fluxMinCut"]
    fluxMaxCut = inputDict["fluxMaxCut"]
    npixelsMinCut = inputDict["npixelsMinCut"]
    npixelsMaxCut = inputDict["npixelsMaxCut"]
    nDonutWanted = inputDict["nDonutWanted"]//2    # divide by 2 for 2 amplifiers
    ellipCut = inputDict["ellipCut"]

    # info about CCDs
    dinfo = decaminfo()

    # time it
    startingTime = time.time()

    
    # list for output of donuts
    donutList = []

    # convert from 2048 by 1024 to 128 by 64, by python magic
    dataBlockView = block_view(dataAmp,(nBlock,nBlock))
    dataBlock = dataBlockView.sum(2).sum(2)

    # find sky background and subtract
    skyBkg = calcSky(dataBlock)
    dataBlock = dataBlock - skyBkg

    # zero edge values - they are screwy
    dataBlock[0,:] = 0.0
    dataBlock[-1,:] = 0.0
    dataBlock[:,-1] = 0.0
    dataBlock[:,0] = 0.0

    # zero pixels lower than the threshold
    dataPeak = numpy.where(dataBlock>fluxThreshold,dataBlock,0.)

    # keep track of which pixels are above threshold too
    # = 1 is > threshold, =0 is below
    # later we'll use this array to keep track of which pixels
    # aren't in a cluster yet
    dataOn = numpy.where(dataBlock>fluxThreshold,1,0)
    
    # find list of pixels that are non-zero
    nzList = numpy.argwhere(dataPeak).tolist()
    numOverThreshold = len(nzList)

    # build a kdTree to locate nearest neighbors
    kdtree = cKDTree(nzList)


    # also sort the array, the [::-1] is python magic to reverse the order of the array!
    ny,nx = dataPeak.shape
    dataPeakFlat = dataPeak.reshape(ny*nx)
    indSort = numpy.argsort(dataPeakFlat)[::-1]

    # get the sorted indices as a tuple of iy,ix 'es
    ixSort = numpy.mod(indSort,nx).tolist()
    iySort = ((indSort - ixSort)//nx).tolist()
    iyxSort = list(zip(iySort,ixSort))

    # now loop over pixels in sorted order
    # but only look at those above the threshold
    for i in range(numOverThreshold):

        # get the iy,ix tuple
        iyx = iyxSort[i]

        # check that this pixel is still available
        if dataOn[iyx]==1:
            
            # try this pixel as the center of a donut
            # find all nearby overthreshold pixels
            d,ind = kdtree.query(iyx,k=kNNFind,distance_upper_bound=distanceLimit)

            # get the centroid of these guys, sum the flux
            # looks like we have to loop
            tempX = numpy.zeros(kNNFind)
            tempY = numpy.zeros(kNNFind)
            tempV = numpy.zeros(kNNFind)
            nOk = 0
            nLengthInd = len(ind)
            for i in range(nLengthInd):            
                # if we don't find enough neighbors, the others have d==inf and ind=nLength
                if ind[i]<numOverThreshold:
                    iy,ix = nzList[ind[i]]
                    # check that this pixel is still available
                    if dataOn[iy,ix] == 1:
                        nOk = nOk + 1
                        tempY[i] = iy
                        tempX[i] = ix
                        tempV[i] = dataPeak[iy,ix]
                        # remove from dataOn array, so we don't reuse
                        dataOn[iy,ix] = 0

            # now all the neighbors are collected
            # calculate flux and centroid and ellipticity, see if its a good guy!
            if nOk>npixelsMinCut:
                # calculate flux, centroid and moments
                flux = tempV.sum()
                xCentroid = (tempV*tempX).sum()/flux
                yCentroid = (tempV*tempY).sum()/flux
                xDiff = tempX - xCentroid
                yDiff = tempY - yCentroid
                xxMoment = (tempV*xDiff*xDiff).sum()/flux
                yyMoment = (tempV*yDiff*yDiff).sum()/flux
                xyMoment = (tempV*xDiff*yDiff).sum()/flux
                if (xxMoment+yyMoment) > 1.e-10:
                    ellip1 = (xxMoment - yyMoment)/(xxMoment + yyMoment)
                    ellip2 = 2.0*xyMoment/(xxMoment+yyMoment)
                else:
                    ellip1 = 0.0
                    ellip2 = 0.0

                # check if this guy passes the cuts

                # calculate the pixels
                ixDonut = int(xCentroid * nBlock) + nBlock//2 + xOffset
                iyDonut = int(yCentroid * nBlock) + nBlock//2 

                # calculate rdecam
                xdecam,ydecam = dinfo.getPosition(extName,ixDonut,iyDonut)
                rdecam = numpy.sqrt(xdecam*xdecam + ydecam*ydecam)
                    
                if nOk<npixelsMaxCut and flux>fluxMinCut and flux<fluxMaxCut and numpy.abs(ellip1)<ellipCut and numpy.abs(ellip2)<ellipCut and rdecam<225.0:
                    # make output list of Dictionaries
                    # convert x,y to 2048x2048 coordinates
                    info = {"x":ixDonut,"y":iyDonut,"mag":flux}
                    donutList.append(info)
                    if len(donutList) >= nDonutWanted:
                        return donutList

                # clear all guys nearby too
#                for j in range(len(ind)):
#                    if d[j]<3.0*distanceLimit:
#                        iy,ix = nzList[ind[j]]
#                        dataOn[iy,ix] = 0

    # Done!
    return donutList



def dfdFinder(hdulist,extName,inputDict={}):
    """  Damn Fast Donut Finder: Find donuts in a F&A CCD.
    The algorithm coalesces the by a factor of N=16 by default.
    Then sorts the list of super-pixels, and loops through each pixel
    above a threshold.  Next it looks for nearest-neighbors inside a radius,
    if some but not too many are present, it sums these, finds their centroid
    and puts them into the output list.
    """

    # cuts and thresholds etc...
    defaultDict = {"nBlock" : 16,
                   "fluxThreshold" : 5.0e3,
                   "kNNFind" : 40,
                   "distanceLimit" : 2.8,
                   "fluxMinCut" : 1.0e4,
                   "fluxMaxCut" : 1.0e8,
                   "npixelsMinCut" : 3,
                   "npixelsMaxCut" : 20,
                   "nDonutWanted" : 20,
                   "ellipCut" : 0.2}

    defaultDict.update(inputDict)


    # grab the data,header
    data = hdulist[extName].data
    header = hdulist[extName].header   

    # time it
    startingTime = time.time()

    # we only want the image portion, drop the overscan
    # HARDCODED!!!
    # and get the two amplifiers separately, since we do no overscan,flat or bias correction
    if extName == "FS4" or extName == "FN4":
        dataAmpB = data[0:2048,56:1080]
        dataAmpA = data[0:2048,1080:2104]
    else:
        dataAmpB = data[50:2098,56:1080]
        dataAmpA = data[50:2098,1080:2104]

    # find em
    donutListA = findDonuts(dataAmpA,1080,defaultDict,extName)
    donutListB = findDonuts(dataAmpB,56,defaultDict,extName)

    # done
    print("dfdFinder took ",time.time()-startingTime," seconds")

    donutList = donutListA + donutListB
    return donutList
    

    
def listToregion(dList,regionFileName="temp.reg"):

    # open region file    
    regf = open(regionFileName,'w')
    regf.write("# Region file format: DS9 version 4.1\n")
    regf.write("# Filename: \n")
    regf.write("global color=green dashlist=8 3 width=1 select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1\n")
    regf.write("\n")

    for donutDict in dList:
        try:
            ix = donutDict["x"]
            iy = donutDict["y"] + 51
            regf.write("image; ellipse(%d,%d,%.5f,%.5f,%.5f)\n" % (ix,iy,10.0,10.0,0.))
        except:
            print(donutDict)

    regf.close()
