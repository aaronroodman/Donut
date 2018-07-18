#
# Utility routines for donut code
#
#     Aaron Roodman (C) SLAC National Accelerator Laboratory, Stanford University 2012.
#
import numpy
import scipy
import scipy.special
import scipy.fftpack as theFFT
import numpy.lib.index_tricks as itricks
import re
import os
from astropy.io import fits as pyfits
import pdb

#
# declare the functions in this file
#
__all__ = ["pupiltopsf","makeXPsf","fouriertrans","invfouriertrans","makeAiry","makeCirc","makeAtmosphere","makeGaussian","makePupilArrays","makeCircWfm","getZemaxArray","getZemaxWfm","getZemaxPsf","fbinrange","flipZemaxArray","loadImage","loadImageFromFile","clipPostageStamp","writePostageStamp","anaPsf","calcStarting","getAllFitsFiles","resample","getFitsWfm"]


#
# Fourier Transform methods
#
def pupiltopsf(xaxis,yaxis,lambdaz,pupil):
#    pupil.resize(xaxis.shape)
# use inverse transform here - needed to get WFM->PSF to match Zemax!!!
    pupilfft = theFFT.ifft2(pupil)
    pupilfftshift = theFFT.fftshift(pupilfft)
    psf = numpy.real(pupilfftshift*numpy.conj(pupilfftshift))
    xpsf,ypsf = makeXPsf(xaxis,yaxis,lambdaz)
    return xpsf,ypsf,psf

def makeXPsf(xaxis,yaxis,lambdaz):

    # see DES logbook Vol1, page 25
    # and page 87

    # allow 1d or 2d arrays
    if xaxis.shape.__len__()==2 :
        xax = xaxis[0,:]
        yax = yaxis[:,0]
    else:
        xax = xaxis
        yax = yaxis

    deltax = xax[1] - xax[0]
    deltay = yax[1] - yax[0]
    Nx = xax.size
    Ny = yax.size
    # now set xpsf and ypsf so that imagearray[Y,X]
    ypsf,xpsf = itricks.mgrid[-1./(2.*deltay):1./(2.*deltay):1./(Ny*deltay),-1./(2.*deltax):1./(2.*deltax):1./(Nx*deltax)]
    xpsf = lambdaz*xpsf
    ypsf = lambdaz*ypsf

    return xpsf,ypsf


def fouriertrans(xaxis,yaxis,input):
    input.resize(xaxis.shape)
    output = theFFT.fft2(input)
    outputshift = theFFT.fftshift(output)
    #
    # see DES logbook Vol1, page 25
    #
    xax = xaxis[0,:]
    yax = yaxis[:,0]

    deltax = xax[1] - xax[0]
    deltay = yax[1] - yax[0]
    Nx = xax.size
    Ny = yax.size
    yout,xout = itricks.mgrid[-1./(2.*deltay):1./(2.*deltay):1./(Ny*deltay),-1./(2.*deltax):1./(2.*deltax):1./(Nx*deltax)]

    return xout,yout,outputshift

def invfouriertrans(xaxis,yaxis,input):
    input.resize(xaxis.shape)
    output = theFFT.ifft2(input)
    outputshift = theFFT.fftshift(output)
    #
    # see DES logbook Vol1, page 25
    #
    xax = xaxis[0,:]
    yax = yaxis[:,0]

    deltax = xax[1] - xax[0]
    deltay = yax[1] - yax[0]
    Nx = xax.size
    Ny = yax.size
    yout,xout = itricks.mgrid[-1./(2.*deltay):1./(2.*deltay):1./(Ny*deltay),-1./(2.*deltax):1./(2.*deltax):1./(Nx*deltax)]

    return xout,yout,outputshift

#
# Make some distributions
#
# use itricks.mgrid[lo:hi:delta] so that axis values go from lo to hi-delta inclusive
# and n=(hi-lo)/delta


def makeAiry(nbin,lo,hi):
    delta = (hi-lo)/nbin
    yaxis,xaxis = itricks.mgrid[lo:hi:delta,lo:hi:delta]
    airyarr = scipy.special.airy(sqrt(xaxis*xaxis+yaxis*yaxis))
    return xaxis,yaxis,airyarr

def makeCirc(nbin,lo,hi,radius):
    delta = (hi-lo)/nbin
    yaxis,xaxis = itricks.mgrid[lo:hi:delta,lo:hi:delta]
    circarray = numpy.where(numpy.sqrt(xaxis*xaxis+yaxis*yaxis)<radius,1.,0.)
    return xaxis,yaxis,circarray

def makeAtmosphere(nbin,binsize,rzero,wavelength,fLength):
    # create the radius array for the Atmosphere
    Lf = 1./binsize
    delta = Lf/nbin
    xAtmos,yAtmos = itricks.mgrid[-Lf/2.:Lf/2.:delta,-Lf/2.:Lf/2.:delta]
    rAtmos = numpy.sqrt(xAtmos*xAtmos + yAtmos*yAtmos)
    arrAtmos = numpy.exp(-3.44*pow(rAtmos*wavelength*fLength/rzero,5./3.))
    xpsfAtmos,ypsfAtmos,psfAtmosC = invfouriertrans(xAtmos,yAtmos,arrAtmos)
    psfAtmos = numpy.abs(psfAtmosC)
    return xpsfAtmos,ypsfAtmos,psfAtmos

def makeGaussian(nbin,lo,hi,sigma):
    delta = (hi-lo)/nbin
    yaxis,xaxis = itricks.mgrid[lo:hi:delta,lo:hi:delta]
    radius2arr = xaxis*xaxis+yaxis*yaxis
    array = numpy.exp(-radius2arr/(2.*sigma*sigma))
    return xaxis,yaxis,array

def makePupilArrays(nbin,lo,hi,radiusOuter):
    # make the rho,theta,xaxis,yaxis,pupilwave,pupilfunc arrays here

    # build the x,y grid, now defined so that imagearray[Y,X]
    delta = (hi-lo)/nbin
    yaxis,xaxis = itricks.mgrid[lo:hi:delta,lo:hi:delta]

    # build the rho,theta areas on the pupil
    rho = numpy.sqrt(xaxis*xaxis+yaxis*yaxis)/radiusOuter
    theta = numpy.arctan2(yaxis,xaxis)

    # pre-create the pupil function array
    dims = rho.shape
    pupilwave = numpy.zeros( dims )     #pupil wavefront phase
    pupilfunc = numpy.zeros( dims ) + 1j * numpy.zeros( dims ) #pupil function

    return xaxis,yaxis,rho,theta,pupilwave,pupilfunc

def makeCircWfm(wfm):
    circarray = numpy.where(wfm!=0.0,1.,0.)
    pupil = circarray*numpy.exp(2j*numpy.pi*wfm)
    return pupil

#
# utilities for Zemax output
#

def getZemaxArray(filename,skiplines,halfsize,nbins):
    #
    # zemax arrays start at xlo,ylo and each line is one value of
    # Y increasing in X, but the first line is the top in Y, and
    # successive lines decrease in Y
    #
    # numpy.loadtxt reads such files into arrays with [Y,X]
    # ordering, so the only flip needed is in Y
    #
    origarray = numpy.loadtxt(filename,skiprows=skiplines)
    fixedarray = numpy.flipud(origarray)

    lo = -halfsize
    hi = halfsize
    delta = halfsize*2.0/nbins
#    xaxis = fbinrange(lo,hi,nbins)
#    yaxis = fbinrange(lo,hi,nbins)

    # for even number of bins, 0.0 is between pixels nbins/2-1 and nbins/2, counting from 0 to nbins-1
    low = lo + delta/2.0
    high = hi - delta/2.0
    yaxis,xaxis = itricks.mgrid[low:high:nbins*1j,low:high:nbins*1j]
    # now xaxis,yaxis have the pixel centers in microns, if halfsize is also in microns
# (Zemax PSF files have their center in the middle of one bin -- as listed in the header
# I believe the center point numbering (ie. row 257 and column 256) is starting from 1)
# right now I am not doing this quite correctly - can fix soon
#


    return xaxis,yaxis,fixedarray

def getZemaxWfm(txtFile,encoding=""):

    # parse for width and number of bins
    regexp1 = re.compile("Exit Pupil Diameter:")
    regexp2 = re.compile("Pupil grid size:")
    regdict = {"width":regexp1,"nbin":regexp2}
    valdict = {"width":0.,"nbin":0}
    zFile = open(txtFile,'r')
    for line in zFile.readlines():
        for varname in list(regdict.keys()):
            aregexp = regdict[varname]
            matobj = aregexp.search(line)
            if matobj != None :
                words = line.split()
                stvalue = words[3]
                value = float(stvalue)
                valdict[varname] = value
                del regdict[varname]
        if len(regdict)==0:
            break
    zFile.close()
    # use meters as the distance unit
    width = valdict["width"] * 1.0e-3
    nbins = int(valdict["nbin"])

    # read in text file
    xaxis,yaxis,data = getZemaxArray(txtFile,16,width/2.0,nbins)

    return xaxis,yaxis,data


def getZemaxPsf(txtFile):

    # parse for width and number of bins
    regexp1 = re.compile("Data area is")
    regexp2 = re.compile("Image grid size:")
    regdict = {"width":regexp1,"nbin":regexp2}
    valdict = {"width":0.,"nbin":0}
    zFile = open(txtFile,'r')
    for line in zFile.readlines():
        for varname in list(regdict.keys()):
            aregexp = regdict[varname]
            matobj = aregexp.search(line)
            if matobj != None :
                words = line.split()
                stvalue = words[3]
                value = float(stvalue)
                valdict[varname] = value
                del regdict[varname]
        if len(regdict)==0:
            break
    zFile.close()
    # use meters as unit
    width = valdict["width"] * 1.0e-6
    nbins = int(valdict["nbin"])

    # read in text file
    xaxis,yaxis,data = getZemaxArray(txtFile,18,width/2.0,nbins)

    return xaxis,yaxis,data


def fbinrange(lo,hi,n):
    binsize = (hi-lo)/n
    locenter = lo + binsize/2.0
    hicenter = hi - binsize/2.0
    vals = numpy.linspace(locenter,hicenter,n)
    return vals


def flipZemaxArray(inarr):
    # reorient Zemax arrays - transpose x and y, and then flip y
    tranarr = inarr.transpose()
    fliparr = numpy.fliplr(tranarr)
    return fliparr



def anaPsf(xaxis,yaxis,data):

    # centroids
    xCentroid = (xaxis*data).sum()/data.sum()
    yCentroid = (yaxis*data).sum()/data.sum()

    # radial sum -- around Centroid
    raxis = numpy.sqrt((xaxis-xCentroid)*(xaxis-xCentroid)+(yaxis-yCentroid)*(yaxis-yCentroid))

    # histogram
    nsumbin = 1000
    npix,bin_edges = scipy.histogram(raxis,nsumbin,(0.,100.))
    rsumpix,bin_edges = scipy.histogram(raxis,nsumbin,(0.,100.),weights=data)

    # calculate ee80
    rsumpixNorm = rsumpix/rsumpix.sum()
    rcumsum = numpy.cumsum(rsumpixNorm)
    # at this point rcumsum[0] is = rsumpixNorm[0], rsumsum[1] to the sum of rsumpixnorm[0] and [1] etc.
    # so rcumsum[0] is the integral for r<bin_edges[1]   and in general rcumsum[i] is the integral of r<bin_edges[i+1]
    # thus icumsum gives the appropriate limits for each rcumsum bin
    #

    icumsum = bin_edges[1:nsumbin+1]
    ee80 = numpy.interp(0.8,rcumsum,icumsum)

    # calculate polarization (w/o seeing, CCD diffusion)
    norm = data.sum()
    qxx = ( (xaxis-xCentroid)*(xaxis-xCentroid)*data ).sum() / norm
    qyy = ( (yaxis-yCentroid)*(yaxis-yCentroid)*data ).sum() / norm
    qyx = ( (yaxis-yCentroid)*(xaxis-xCentroid)*data ).sum() / norm

    e1 = (qxx-qyy)/(qxx+qyy)
    e2 = 2.0*qyx/(qxx+qyy)

    return ee80,qxx,qyy,qyx,e1,e2

def loadImageFromFile(fileName,iExtension=0,ixCenter=32,iyCenter=32,nhalfPixels=32,dumpFlag=True,gain=1.0,constantError=7.1,outputPrefix="test"):

    # nPixels
    nPixels = 2*nhalfPixels

    # open fits file
    hdulist = pyfits.open(fileName)
    if dumpFlag:
        hdulist.info()
    data = hdulist[iExtension].data

    # get postage stamp:
    # designed for even dimensioned array, where the center bin is the one just after the central point
    # so if nhalfPixel = 32 and ixCenter=32,iyCenter=32   xlo = 0 yhi = 63
    #
    # note that pyfits arrays are [row,column] or [Y,X]
    xlo = int(ixCenter - nhalfPixels)
    xhi = int(ixCenter + nhalfPixels - 1)
    ylo = int(iyCenter - nhalfPixels)
    yhi = int(iyCenter + nhalfPixels - 1)

    # check that xlo,xhi  ylo,yhi are all within range for data
    # if not then clip the array, but place it, centered, in an array of the desired dimension
    xloLimit = 0
    xhiLimit = data.shape[1] - 1
    yloLimit = 0
    yhiLimit = data.shape[0] - 1

    xlopad = 0
    xhipad = 0
    ylopad = 0
    yhipad = 0

    if xlo<xloLimit:
        xlopad = xloLimit - xlo
        xlo = xloLimit
    if xhi>xhiLimit:
        xhipad = xhi - xhiLimit
        xhi = xhiLimit
    if ylo<yloLimit:
        ylopad = yloLimit - ylo
        ylo = yloLimit
    if yhi>yhiLimit:
        yhipad = yhi - yhiLimit
        yhi = yhiLimit

    # imgarray is in electrons
    # data[ylo:yhi+1] gets bins from ylo to yhi
    if xlopad==0 and xhipad == 0 and ylopad == 0 and yhipad == 0:
        imgarray = gain * data[ylo:yhi+1,xlo:xhi+1]
    else:
        tmparray = gain * data[ylo:yhi+1,xlo:xhi+1]
        imgarray = numpy.zeros((nPixels,nPixels))
        imgarray[ylopad:nPixels-yhipad,xlopad:nPixels-xhipad] = tmparray

    # dump postage stamp
    if dumpFlag:
        hdu = pyfits.PrimaryHDU(imgarray)
        hdulist = pyfits.HDUList([hdu])
        outFile = outputPrefix + "-image.fits"
        hdulist.writeto(outFile,overwrite=True)
        # close file
        hdulist.close()

    weightsq = imgarray + constantError*constantError

    return imgarray,weightsq

def loadImage(hdulist,iExtension=0,ixCenter=32,iyCenter=32,nhalfPixels=32,dumpFlag=True,gain=1.0,constantError=7.1,outputPrefix="test"):

    # nPixels
    nPixels = 2*nhalfPixels

    # get data from fits file
    if dumpFlag:
        hdulist.info()
    data = hdulist[iExtension].data

    # get postage stamp:
    # designed for even dimensioned array, where the center bin is the one just after the central point
    # so if nhalfPixel = 32 and ixCenter=32,iyCenter=32   xlo = 0 yhi = 63
    #
    # note that pyfits arrays are [row,column] or [Y,X]
    xlo = int(ixCenter - nhalfPixels)
    xhi = int(ixCenter + nhalfPixels - 1)
    ylo = int(iyCenter - nhalfPixels)
    yhi = int(iyCenter + nhalfPixels - 1)

    # check that xlo,xhi  ylo,yhi are all within range for data
    # if not then clip the array, but place it, centered, in an array of the desired dimension
    xloLimit = 0
    xhiLimit = data.shape[1] - 1
    yloLimit = 0
    yhiLimit = data.shape[0] - 1

    xlopad = 0
    xhipad = 0
    ylopad = 0
    yhipad = 0
    if xlo<xloLimit:
        xlopad = xloLimit - xlo
        xlo = xloLimit
    if xhi>xhiLimit:
        xhipad = xhi - xhiLimit
        xhi = xhiLimit
    if ylo<yloLimit:
        ylopad = yloLimit - ylo
        ylo = yloLimit
    if yhi>yhiLimit:
        yhipad = yhi - yhiLimit
        yhi = yhiLimit

    # imgarray is in electrons
    # data[ylo:yhi+1] gets bins from ylo to yhi
    if xlopad==0 and xhipad == 0 and ylopad == 0 and yhipad == 0:
        imgarray = gain * data[ylo:yhi+1,xlo:xhi+1]
    else:
        tmparray = gain * data[ylo:yhi+1,xlo:xhi+1]
        imgarray = numpy.zeros((nPixels,nPixels))
        imgarray[ylopad:nPixels-yhipad,xlopad:nPixels-xhipad] = tmparray

    # dump postage stamp
    if dumpFlag:
        hdu = pyfits.PrimaryHDU(imgarray)
        hdulistOut = pyfits.HDUList([hdu])
        outFile = outputPrefix + "-image.fits"
        hdulistOut.writeto(outFile,overwrite=True)
        # close file
        hdulistOut.close()

    weightsq = imgarray + constantError*constantError

    return imgarray,weightsq


def clipPostageStamp(ixCenter,iyCenter,data,nPixels):
    # clips out a subarray - and works even if the postage stamp runs off the edge

    nhalfPixels = nPixels/2
    # designed for even dimensioned array, where the center bin is the one
    #               just after the central point
    # so if nhalfPixel = 32 and ixCenter=32,iyCenter=32   xlo = 0 yhi = 63
    #
    # note that pyfits arrays are [row,column] or [Y,X]
    xlo = int(ixCenter - nhalfPixels)
    xhi = int(ixCenter + nhalfPixels - 1)
    ylo = int(iyCenter - nhalfPixels)
    yhi = int(iyCenter + nhalfPixels - 1)

    # check that xlo,xhi  ylo,yhi are all within range for data
    # if not then clip the array, but place it, centered, in an array of the desired dimension
    xloLimit = 0
    xhiLimit = data.shape[1] - 1
    yloLimit = 0
    yhiLimit = data.shape[0] - 1

    xlopad = 0
    xhipad = 0
    ylopad = 0
    yhipad = 0
    if xlo<xloLimit:
        xlopad = xloLimit - xlo
        xlo = xloLimit
    if xhi>xhiLimit:
        xhipad = xhi - xhiLimit
        xhi = xhiLimit
    if ylo<yloLimit:
        ylopad = yloLimit - ylo
        ylo = yloLimit
    if yhi>yhiLimit:
        yhipad = yhi - yhiLimit
        yhi = yhiLimit

    # imgarray is in electrons
    # data[ylo:yhi+1] gets bins from ylo to yhi
    if xlopad==0 and xhipad == 0 and ylopad == 0 and yhipad == 0:
        stamp = data[ylo:yhi+1,xlo:xhi+1]
    else:
        tmparray = data[ylo:yhi+1,xlo:xhi+1]
        stamp = numpy.zeros((nPixels,nPixels))
        stamp[ylopad:nPixels-yhipad,xlopad:nPixels-xhipad] = tmparray

    # return the stamp
    return stamp



def writePostageStamp(hdulist,iFile,iExtension=0,ixCenter=32,iyCenter=32,nhalfPixels=32,outputFile="temp",extraheader=None):

    # nPixels
    nPixels = 2*nhalfPixels

    # get pointer to this Extension's data
    data = hdulist[iExtension].data

    # get postage stamp:
    imgarray = clipPostageStamp(ixCenter,iyCenter,data,nhalfPixels)

    # dump postage stamp to new file
    stamphdu = pyfits.PrimaryHDU(imgarray)
    stamphdulist = pyfits.HDUList([stamphdu])

    # fill some header words of the postage stamp
    stamphdr = stamphdu.header
    stamphdr.update("IFILE",iFile ,"File Number")
    stamphdr.update("IEXT",iExtension,"Fits file extension")
    stamphdr.update("IX",int(ixCenter),"nominal x Center pixel")
    stamphdr.update("IY",int(iyCenter),"nominal y Center pixel")

    # add extra header words
    if extraheader!=None:
        for key in extraheader:
            stamphdr.update(key,extraheader[key])

    # write out postage stamp
    stamphdulist.writeto(outputFile,overwrite=True)


def calcWeight(imgarray,constantError):
    weightsq = (constantError*constantError)/(4.0*imgarray) + 0.25
    weight = numpy.sqrt(weightsq)



def calcStarting(imgarray,nsigma=1.5,printLevel=0,debugFlag=False):

    # copy image to temporary array
    locarray = imgarray.copy()

    # size of postage stamp
    nArea = locarray.shape[0]*locarray.shape[1]

    # get Mean,Sigma, Mask of pixels below Mean+nsigma*sigma
    # use numpy Masked Array!
    imgMean = locarray.mean()
    imgStd = locarray.std()
    if debugFlag:
        print("calcStarting: mean,sigma = ",imgMean,imgStd)
    imgarrayMask = numpy.ma.masked_greater(locarray,imgMean+1*imgStd)
    countsNew = (imgarrayMask.mask==False).sum()
    countsOld = -1

    while countsOld!=countsNew:
        countsOld = countsNew
        imgMean = imgarrayMask.mean()
        imgStd = imgarrayMask.std()
        if printLevel>=2:
            print("calcStarting: mean,sigma,counts = ",imgMean,imgStd,countsOld)
        imgarrayMask = numpy.ma.masked_greater(locarray,imgMean+nsigma*imgStd)
        countsNew = (imgarrayMask.mask==False).sum()

        ## # now find x,y-local mean of pixels above pedestal
        ## yloc,xloc = itricks.mgrid[-nPixBig:nPixBig+1,-nPixBig:nPixBig+1]
        ## locarraySubt = locarray - imgMean
        ## total = locarraySubt.sum()
        ## xLocMean = (locarraySubt*xloc).sum()/total
        ## yLocMean = (locarraySubt*yloc).sum()/total

        ## # next adjust the iX and iY to be centered on the local position
        ## ixNewCenter = ixCenter + int(xLocMean)
        ## iyNewCenter = iyCenter + int(yLocMean)
        ## self.loadImage(filename,iExtension,ixNewCenter,iyNewCenter,self.nPixels,True)

        # return Pedestal value and Sum of values above the pedestal
    return imgMean,locarray.sum()-nArea*imgMean


def getAllFitsFiles(path):
    """ given a directory, this method finds all fits files and returns lists with their name and number
    """
    names = []
    numbers = []

    # list all files
    files = os.listdir(path)

    # ASSUME name is of form TEXTNUMBER
    pattern = "(\D+)(\d+)"
    prog = re.compile(pattern)

    # loop over files
    for file in files:
        basename = os.path.basename(file)
        base, exten = os.path.splitext(basename)
        # only use .fits files
        if exten == ".fits":
            names.append(basename)
            # find number
            result = prog.match(base)
            text,num = result.groups()
            numbers.append(int(num))

    # all done
    return names,numbers

def resample(inArr,inPixscale,outShape,outPixscale):
    """ resamples the inArr with inPixscale to outShape and outPixscale
    this method does NO interpolation, just samples the input array
    """

    innRow,innCol = inArr.shape
    outnRow,outnCol = outShape
    
    # construct map from outputRow to inputRow
    iRowBinEdges = numpy.zeros((innRow+1))   # bin edges 
    oRowBinCenters = numpy.zeros((outnRow))

    iRowCenter = (innRow-1.)/2.   #center in pixel numbering, ie. center=511.5 for n=1024 and i=0->0123
    for iRow in range(innRow+1):
        iRowBinEdges[iRow] = (iRow - iRowCenter - 0.5)*inPixscale

    oRowCenter = (outnRow-1.)/2.   #center in pixel numbering, ie. center=511.5 for n=1024 and i=0->0123
    for oRow in range(outnRow):
        oRowBinCenters[oRow] = (oRow - oRowCenter)*outPixscale

    # find output centers in input bins, 
    oRowMap = numpy.digitize(oRowBinCenters,iRowBinEdges) - 1

    # construct map from outputCol to inputCol
    iColBinEdges = numpy.zeros((innCol+1))   # bin edges 
    oColBinCenters = numpy.zeros((outnCol))

    iColCenter = (innCol-1.)/2.   #center in pixel numbering, ie. center=511.5 for n=1024 and i=0->0123
    for iCol in range(innCol+1):
        iColBinEdges[iCol] = (iCol - iColCenter - 0.5)*inPixscale

    oColCenter = (outnCol-1.)/2.   #center in pixel numbering, ie. center=511.5 for n=1024 and i=0->0123
    for oCol in range(outnCol):
        oColBinCenters[oCol] = (oCol - oColCenter)*outPixscale

    # find output centers in input bins
    oColMap = numpy.digitize(oColBinCenters,iColBinEdges) - 1

    # now resample...
    outArr = numpy.zeros((outnRow,outnCol))
    for orow in range(outnRow):
        for ocol in range(outnCol):
            irow = oRowMap[orow]
            icol = oColMap[ocol]
            if irow>=0 and irow<innRow and  icol>=0 and icol<innCol:
                outArr[orow,ocol] = inArr[irow,icol]

    # done
    return outArr


def getFitsWfm(fitsFile,extNo):

    # parse for width and number of bins
    hdu = pyfits.open(fitsFile)
    data = hdu[extNo].data
    data64 = data.astype(numpy.float64)
    return data64
    
    
