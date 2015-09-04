import numpy
import scipy
import time
import pyfits
import copy
import statsmodels.api as sm
import os.path
from donutlib.PointMesh import PointMesh
from decamutil import decaminfo 
from decamutil import mosaicinfo
try:
    from ROOT import TTree, TFile, gROOT, TCanvas, gStyle, TGraph2D, TGraphErrors, SetOwnership
    from hplotlib import hfillhist
except:
    print "donutana: could not import ROOT"
import pdb

class donutana(object):
    """ donutana is a class used to analyze the results from donut fits for the DES experiment.

    Aaron Roodman (C) SLAC National Accelerator Laboratory, Stanford University 2012.
    """

    def __init__(self,**inputDict):

        # default values for calibration information

        #  see doit.txt file 10/14/2012
        #  comes from 20120914seq1 suplmented by 20121001seq1 for tilts
        #  I've rounded to 2 sig. figures, removed elements to correct for 15deg rotation,
        #  symmetrized, then added the LR sub-matrix to match Zemax,
        #  then finally inverted and rounded.
        #  
        #  hexapodM = 
        #      matrix ([[  0.0e+00,  -2.6e+05,   4.5e+03,   0.0e+00],
	#               [ -2.6e+05,  -0.0e+00,  -0.0e+00,  -4.5e+03],
       	#               [ -5.3e+04,  -0.0e+00,  -0.0e+00,  -9.8e+01],
        #               [  0.0e+00,   5.3e+04,  -9.8e+01,   0.0e+00]])

        hexapodArrayDiagonal = numpy.array(((  0.0e+00,  -2.6e+05,   4.5e+03,   0.0e+00),
                                    ( -2.6e+05,   0.0e+00,   0.0e+00,  -4.5e+03),
                                    ( -5.3e+04,   0.0e+00,   0.0e+00,  -9.8e+01),
                                    (  0.0e+00,   5.3e+04,  -9.8e+01,   0.0e+00)))

        # now replace this with the array from the actual measurements, which
        # still includes the 15deg hexapod rotation

        hexapodArrayRotated =  numpy.array(((-1.9e5, -1.6e5,  4.1e3, -1.3e3),
                                     (-1.6e5,  1.9e5, -1.3e3, -4.1e3),
                                     (-4.8e4,  1.6e4,  -30. ,  -88.),
                                     ( 1.6e4,  4.8e4,  -88. ,  30. )))

        hexapodArray20121020 = numpy.array(((  0.00e+00 ,  1.07e+05 ,  4.54e+03 ,  0.00e+00),
                                            (  1.18e+05 , -0.00e+00 ,  0.00e+00 , -4.20e+03),
                                            ( -4.36e+04 ,  0.00e+00 ,  0.00e+00 , -8.20e+01),
                                            (  0.00e+00 ,  4.42e+04 , -8.10e+01 ,  0.00e+00) ))
                 
        self.paramDict  = {"z4Conversion":172.0,
                           "deltaZSetPoint":0.0,
                           "alignmentMatrix":hexapodArray20121020,
                           "zPointsFile":"",
                           "z4PointsFile":"",
                           "z5PointsFile":"",
                           "z6PointsFile":"",
                           "z7PointsFile":"",
                           "z8PointsFile":"",
                           "z9PointsFile":"",
                           "z10PointsFile":"",
                           "z11PointsFile":"",
                           "z12PointsFile":"",
                           "z13PointsFile":"",
                           "z14PointsFile":"",
                           "z15PointsFile":"",
                           "nInterpGrid":8,
                           "interpMethod":"idw",
                           "methodVal":None,
                           "donutCutString":"",
                           "sensorSet":"FandAOnly",    # or "ScienceOnly" or "Mosaic"
                           "unVignettedOnly":False,
                           "doTrefoil":False,
                           "doSpherical":False,
                           "doQuadrefoil":False,
                           "doRzero":False,
                           "histFlag":False,
                           "debugFlag":False}
        
        # update paramDict from inputDict
        self.paramDict.update(inputDict)

        # get DECam geometry information
        if self.paramDict["sensorSet"] == "ScienceOnly" or self.paramDict["sensorSet"] == "FandAOnly" or self.paramDict["sensorSet"] == "PlusOnly" or self.paramDict["sensorSet"] == "MinusOnly"  :
            self.infoObj = decaminfo()
        else:
            print "HEY FOOL, set either ScienceOnly or FandAOnly !!!!"
            exit()
            self.infoObj = mosaicinfo()
            
        self.info = self.infoObj.infoDict

        # fill PointMesh coordinate list and gridDict
        self.coordList = []
        self.gridDict = {}

        #  Code for crude vignetting cut - now obsolete
        if self.paramDict["unVignettedOnly"]:
            limitsFA = {"FS1":[0,1024],"FS2":[0,1024],"FS3":[0,1024],"FS4":[0,1024],
                        "FN1":[-1024,0],"FN2":[-1024,0],"FN3":[-1024,0],"FN4":[-1024,0]}
        else:
            limitsFA = {"FS1":[-1024,1024],"FS2":[-1024,1024],"FS3":[-1024,1024],"FS4":[-1024,1024],
                        "FN1":[-1024,1024],"FN2":[-1024,1024],"FN3":[-1024,1024],"FN4":[-1024,1024]}

        # build list of ccds - options are FandAOnly, ScienceOnly, All
        if self.paramDict["sensorSet"] == "FandAOnly" :
            for ccd in self.info.keys():
                ccdinfo = self.info[ccd]
                if  ccdinfo["FAflag"]:
                    self.coordList.append(ccd)
        elif self.paramDict["sensorSet"] == "PlusOnly" :
            for ccd in self.info.keys():
                ccdinfo = self.info[ccd]
                if  ccdinfo["FAflag"] and ccdinfo["Offset"] == 1500.0:
                    self.coordList.append(ccd)
        elif self.paramDict["sensorSet"] == "MinusOnly" :
            for ccd in self.info.keys():
                ccdinfo = self.info[ccd]
                if  ccdinfo["FAflag"] and ccdinfo["Offset"] == -1500.0:
                    self.coordList.append(ccd)
        elif self.paramDict["sensorSet"] == "ScienceOnly":
            for ccd in self.info.keys():
                ccdinfo = self.info[ccd]
                if  not ccdinfo["FAflag"]:
                    self.coordList.append(ccd)
        elif  self.paramDict["sensorSet"] == "Both":
            for ccd in self.info.keys():
                self.coordList.append(ccd)
        else:
            print "donutana.py: HEY FOOL, you have to specify FandAOnly, ScienceOnly, or Both"
            exit()

        # reset doRzero based on sensorSet - why do this?
        #if self.paramDict["sensorSet"] != "FandAOnly":
        #    self.paramDict["doRzero"] = False

        # loop over ccds
        for ccd in self.coordList:
            ccdinfo = self.info[ccd]
            # either FA sensors or Science sensors
            if  ( ccdinfo["FAflag"] ):

                xlo = ccdinfo["xCenter"] - 1024 * self.infoObj.mmperpixel
                xhi = ccdinfo["xCenter"] + 1024 * self.infoObj.mmperpixel
                
                ylo = ccdinfo["yCenter"] - 1024 * self.infoObj.mmperpixel
                yhi = ccdinfo["yCenter"] + 1024 * self.infoObj.mmperpixel

                # fill gridDict
                self.gridDict[ccd] = [self.paramDict["nInterpGrid"],ylo,yhi,self.paramDict["nInterpGrid"],xlo,xhi]

            elif ( self.paramDict["sensorSet"] == "ScienceOnly" and not ccdinfo["FAflag"] ):

                xlo = ccdinfo["xCenter"] - 1024 * self.infoObj.mmperpixel
                xhi = ccdinfo["xCenter"] + 1024 * self.infoObj.mmperpixel
                ylo = ccdinfo["yCenter"] - 2048 * self.infoObj.mmperpixel
                yhi = ccdinfo["yCenter"] + 2048 * self.infoObj.mmperpixel

                # fill gridDict
                # 1/31/2014 - make the grid x2 bigger in Y than X, to match CCD sizes!
                self.gridDict[ccd] = [2*self.paramDict["nInterpGrid"],ylo,yhi,self.paramDict["nInterpGrid"],xlo,xhi]

            #elif self.paramDict["sensorSet"] == "Mosaic" :
            #    xlo = ccdinfo["xCenter"] - 1024 * self.infoObj.mmperpixel
            #    xhi = ccdinfo["xCenter"] + 1024 * self.infoObj.mmperpixel
            #    ylo = ccdinfo["yCenter"] - 2048 * self.infoObj.mmperpixel
            #    yhi = ccdinfo["yCenter"] + 2048 * self.infoObj.mmperpixel
            #
            #    # fill gridDict
            #    self.gridDict[ccd] = [self.paramDict["nInterpGrid"],ylo,yhi,self.paramDict["nInterpGrid"],xlo,xhi]


        # also keep a dictionary to the meshes
        self.meshDict = {}

        # for backward compatibility...
        if self.paramDict["zPointsFile"] != "":
            self.paramDict["z4PointsFile"] = self.paramDict["zPointsFile"]

        # build the reference meshes
        for iZ in range(4,15+1):
            name = "z%dPointsFile" % (iZ)
            title = "Zernike %d" % (iZ)
            if self.paramDict[name] != "":
                # check that the file exists!
                if os.path.isfile(self.paramDict[name]):
                    theMesh = PointMesh(self.coordList,self.gridDict,pointsFile=self.paramDict[name],myMethod=self.paramDict["interpMethod"],methodVal=self.paramDict["methodVal"],title=title)
                    meshName = "z%dMesh" % iZ
                    self.meshDict[meshName] = theMesh

        # matrix for Hexapod calculation
        self.alignmentMatrix = numpy.matrix(self.paramDict["alignmentMatrix"])
        self.alignmentMatrixSquared = numpy.matrix(self.paramDict["alignmentMatrix"]*self.paramDict["alignmentMatrix"])

        # store some constants
        arcsecperrad = 3600.0 * 180. / numpy.pi
        mmpermicron = 0.001
        self.zangleconv = arcsecperrad * mmpermicron

    def fillPoints(self,dataDict,extraCut=""):

        # this fills a dictionary, keyed by Zernike term (and rzero,chi2,nele), of dictionaries, keyed by Coord, of lists (which are converted to arrays)

        bigDict = {}

        # make Points dictionaries for PointMeshes from the list of Donut dictionaries
        for iZ in range(4,8+1):
            zkey = "z%d" % (iZ)
            bigDict[zkey] = {}
        if self.paramDict["doTrefoil"]:
            bigDict["z9"] = {}
            bigDict["z10"] = {}
        if self.paramDict["doSpherical"]:
            bigDict["z11"] = {}
        if self.paramDict["doQuadrefoil"]:
            bigDict["z14"] = {}
            bigDict["z15"] = {}
        if self.paramDict["doRzero"]:
            bigDict["rzero"] = {}
            bigDict["chi2"] = {}
            bigDict["nele"] = {}

        # make a blank list for each Coord for each of the Pnts dicts
        for key in bigDict.keys():
            for coord in self.coordList:
                bigDict[key][coord] = []
                        
        # fill bigDict lists from Donut information dictionaries
        for donut in dataDict:

            if type(donut) == dict or type(donut) == pyfits.header.Header:

                try:
                    extname = donut["EXTNAME"]
                except:
                    extname = str(donut["IEXT"]-1) 
                    
                ifile = donut["IFILE"]
                zern4 = donut["ZERN4"]
                zern5 = donut["ZERN5"]
                zern6 = donut["ZERN6"]
                zern7 = donut["ZERN7"]
                zern8 = donut["ZERN8"]

                if self.paramDict["doTrefoil"]:
                    zern9 = donut["ZERN9"]
                    zern10 = donut["ZERN10"]

                if self.paramDict["doSpherical"]:
                    zern11 = donut["ZERN11"]

                if self.paramDict["doQuadrefoil"]:
                    zern14 = donut["ZERN14"]
                    zern15 = donut["ZERN15"]

                if self.paramDict["doRzero"]:
                    rzero = donut["rzero"]
                    
                ix = donut["IX"]
                iy = donut["IY"]
                chi2 = donut["CHI2"]
##                fitstat = donut["FITSTAT"]   ## why commented out?
                nele = donut["NELE"]
##                bkgd = donut["BKGD"]
                
            else:
                
#                if type(donut.iext) == float or type(donut.iext) == int :
#                    extname = str(int(donut.iext)-1)   #KLUDGE for MOSAIC
#                elif type(donut.iext) == str:
#                    extname = donut.iext

                extname = donut.extname
                try:
                    ifile = donut.ifile
                except:
                    ifile = 0

                # for focus scans
                try:
                    ifocus = donut.ifocus
                except:
                    ifocus = 0
                    
                zern4 = donut.zern4
                zern5 = donut.zern5 
                zern6 = donut.zern6 
                zern7 = donut.zern7 
                zern8 = donut.zern8 
                zern9 = donut.zern9
                zern10 = donut.zern10
                ix = donut.ix
                iy = donut.iy

                try:
                    chi2 = donut.chi2
                    fitstat = donut.fitstat
                    nele = donut.nele
                    bkgd = donut.bkgd
                except:
                    chi2 = 1.
                    fitstat = 1
                    nele = 1.
                    bkgd = 0.

                try:
                    zern11 = donut.zern11
                except:
                    zern11 = 0.

                try:
                    zern14 = donut.zern14
                    zern15 = donut.zern15
                except:
                    zern14 = 0.
                    zern15 = 0.

                try:
                    rzero = donut.rzero
                except:
                    rzero = 0.

            # compound cut for 2 bad amplifiers needed
            ampOk = True
            if extname=="FS4" and ix>1024 and numpy.log10(nele)>5.7:
                ampOk = False
            if extname=="FN2" and ix<1024 and numpy.log10(nele)>5.9:
                ampOk = False

            # apply cut (could use compile command here)
            good = False
            if self.paramDict["donutCutString"]=="" and extraCut=="":
                good = True
            elif self.paramDict["donutCutString"]!="" and extraCut=="":
                if eval(self.paramDict["donutCutString"]):
                    good = True
            else:
                if eval(self.paramDict["donutCutString"] + " and " + extraCut):
                    good = True

            # good Donut
            if good and bigDict["z4"].has_key(extname):
                xval,yval = self.infoObj.getPosition(extname,ix,iy)
                # convert from z4 to dz
                dz = zern4 * self.paramDict["z4Conversion"]
                wgt = 1.0

                # append to relevant list
                # make dicts of pnts
                pntsDict = {}
                pntsDict["z4"] = [xval,yval,dz,wgt]
                pntsDict["z5"] = [xval,yval,zern5,wgt]
                pntsDict["z6"] = [xval,yval,zern6,wgt]
                pntsDict["z7"] = [xval,yval,zern7,wgt]
                pntsDict["z8"] = [xval,yval,zern8,wgt]
                if self.paramDict["doTrefoil"]:
                    pntsDict["z9"] = [xval,yval,zern9,wgt]
                    pntsDict["z10"] = [xval,yval,zern10,wgt]

                if self.paramDict["doSpherical"]:
                    pntsDict["z11"] = [xval,yval,zern11,wgt]

                if self.paramDict["doQuadrefoil"]:
                    pntsDict["z14"] = [xval,yval,zern14,wgt]
                    pntsDict["z15"] = [xval,yval,zern15,wgt]

                if self.paramDict["doRzero"]:
                    pntsDict["rzero"] = [xval,yval,rzero,wgt]
                    pntsDict["chi2"] = [xval,yval,chi2,wgt]
                    pntsDict["nele"] = [xval,yval,nele,wgt]

                # now enter these points into bigDict
                for key in bigDict.keys():
                    bigDict[key][extname].append(pntsDict[key])
                    


        # convert all the lists to numpy arrays
        allPointsDict = {}
        for key in bigDict.keys():
            allPointsDict[key] = {}
            for coord in self.coordList:
                allPointsDict[key][coord] = numpy.array(bigDict[key][coord])
     
        return allPointsDict
        
    def makeMeshes(self,donutData,extraCut="",myMethod='none',methodVal=None):
        # make the meshes from the data

        # fill dictionaries of Points, separated by extension
        allPointsDict = self.fillPoints(donutData,extraCut)

        # build meshes for all keys present
        allMeshDict = {}
        for key in allPointsDict.keys():

            thePoints = allPointsDict[key]
            theMesh = PointMesh(self.coordList,self.gridDict,pointsArray=thePoints,myMethod=myMethod,methodVal=methodVal)
            meshName = "%sMesh" % (key)
            allMeshDict[meshName] = theMesh

        return allMeshDict

    def analyzeDonuts(self,donutData,extraCut="",doCull=False,cullCut=0.90):
        """ analyzeDonuts takes a list of dictionaries with donut information as input, containing the results of the donut fits
        """
        
        dictOfMeshes = self.makeMeshes(donutData,extraCut=extraCut)
        donutDict = self.analyzeMeshes(dictOfMeshes,doCull=doCull,cullCut=cullCut)
        return donutDict


    def analyzeMeshes(self,dictOfMeshes,doCull=False,cullCut=0.90):
        """ analyzeMeshes takes a dictionary with meshes as input, and fits it to the references
        output dictionary includes a deepcopy of the input meshes
        note:
           this fits to zDiff = dictOfMeshes - interpolation of self.meshDict
           ie. the sign is defined as (other guy) - (me)
        """

        # deepcopy the meshes! no longer adjusts in place  
        dictOfMeshesCopy = {}
        for key in dictOfMeshes.keys():
            theMesh = copy.deepcopy(dictOfMeshes[key])
            dictOfMeshesCopy[key] = theMesh

        # then we analyze the Zernike polynomials by comparing them to stored references
        # analyzeDonuts returns a dictionary containing deltas and rotation angles for focus, astigmatism and coma
        # as well as the hexapod adjustments
        # it also makes a Canvas of plots for study

        dictOfResults = {}

        # loop over keys, fitting References
        for key in dictOfMeshesCopy.keys():

            # some special cases
            resultsKeyName = "%sResultDict" % (key.replace("Mesh",""))
            meshName = key
            if key=="z4Mesh":
                dictOfResults[resultsKeyName] = self.fitToRefMesh(self.meshDict[meshName],dictOfMeshesCopy[meshName],self.zangleconv)
            elif key=="rzeroMesh":
                dictOfResults[resultsKeyName] = self.analyzeRzero(dictOfMeshesCopy["rzeroMesh"],dictOfMeshesCopy["chi2Mesh"],dictOfMeshesCopy["neleMesh"])
            elif key=="chi2Mesh":
                continue
            elif key=="neleMesh":
                continue
            else:
                # check that meshDict has the appropriate mesh
                if self.meshDict.has_key(meshName):
                    dictOfResults[resultsKeyName] = self.fitToRefMesh(self.meshDict[meshName],dictOfMeshesCopy[meshName])

        # if we want to cull, cull and refit!
        # use the fitted weight to cull - culltype="fit" 
        if doCull:
            # this will take the wgt's from the fit and take their product
            # for each point, and cull at a product of cullCut
            self.cullAllMeshes(dictOfMeshesCopy,cullCut=cullCut)

            for key in dictOfMeshesCopy.keys():

                # some special cases
                resultsKeyName = "%sResultDict" % (key.replace("Mesh",""))
                meshName = key
                if key=="z4Mesh":
                    dictOfResults[resultsKeyName] = self.fitToRefMesh(self.meshDict[meshName],dictOfMeshesCopy[meshName],self.zangleconv)
                elif key=="rzeroMesh":
                    dictOfResults[resultsKeyName] = self.analyzeRzero(dictOfMeshesCopy["rzeroMesh"],dictOfMeshesCopy["chi2Mesh"],dictOfMeshesCopy["neleMesh"])
                elif key=="chi2Mesh":
                    continue
                elif key=="neleMesh":
                    continue
                else:
                    # check that meshDict has the appropriate mesh
                    if self.meshDict.has_key(meshName):
                        dictOfResults[resultsKeyName] = self.fitToRefMesh(self.meshDict[meshName],dictOfMeshesCopy[meshName])
            
        # analyze this data and extract the hexapod coefficients
        donutDict = self.calcHexapod(dictOfResults)

        if len(donutDict)==0:
            goodCalc = False
        else:
            goodCalc = True

        # add the individual fit results here too
        for key in dictOfResults.keys():
            donutDict[key] = dictOfResults[key]

        # and add the meshes too
        for key in dictOfMeshesCopy.keys():
            donutDict[key] = dictOfMeshesCopy[key]


        # make a Canvas of plots for this image
        # plot Histogram of Difference before fit, after fit, and after fit vs. X,Y position
        if self.paramDict["histFlag"] and dictOfResults["z4ResultDict"].has_key("deltaArrayBefore") and goodCalc:

            # setup plots
            gStyle.SetStatH(0.32)
            gStyle.SetStatW(0.4)
            gStyle.SetOptStat(1111111)
            gStyle.SetMarkerStyle(20)
            gStyle.SetMarkerSize(0.5)
            gStyle.SetPalette(1)            
            gROOT.ForceStyle()
            
            # loop over results, making plots for each
            nplots = 0
            plotDict = {} 
            for key in dictOfResults.keys():
                theResultDict = dictOfResults[key]

                keyId = key.replace("ResultDict","")
                # special cases 
                if keyId=="z4":
                    nWavesBefore = 200.0
                    nWavesAfter = 200.0
                else:
                    nWavesBefore = 1.0
                    nWavesAfter = 0.2

                if keyId!="rzero":
                    nplots = nplots + 1
                    hBefore = hfillhist(key+"Before","Delta "+key+", Before Fit",theResultDict["deltaArrayBefore"],200,-nWavesBefore,nWavesBefore)
                    hAfter = hfillhist(key+"zAfter","Delta "+key+", After Fit",theResultDict["deltaArrayAfter"],200,-nWavesAfter,nWavesAfter)
                    hBefore2D = TGraph2D(key+"Before2D","Delta "+key+", Before Fit, vs. Position;X[mm];Y[mm]",theResultDict["deltaArrayBefore"].shape[0],theResultDict["deltaArrayX"],theResultDict["deltaArrayY"],theResultDict["deltaArrayBefore"])
                    hAfter2D = TGraph2D(key+"After2D","Delta "+key+", After Fit, vs. Position;X[mm];Y[mm]",theResultDict["deltaArrayAfter"].shape[0],theResultDict["deltaArrayX"],theResultDict["deltaArrayY"],theResultDict["deltaArrayAfter"])
                   
                    plotList = [hBefore,hAfter,hBefore2D,hAfter2D]
                    plotDict[keyId] = plotList

            # the Canvas

            # unique name for our canvas
            tstr = "canvas" + str(time.time())

            canvas = TCanvas(tstr,tstr,300*nplots,1000)
            canvas.Divide(nplots,4)

            # plot em
            jZ = 0
            for iZ in range(4,15+1):
                key = "z%d" % (iZ)
                if plotDict.has_key(key):
                    jZ = jZ + 1

                    plotList = plotDict[key]

                    icanvas = jZ + 0*nplots
                    canvas.cd(icanvas)
                    plotList[0].Draw()
                    icanvas = jZ + 1*nplots
                    canvas.cd(icanvas)
                    plotList[1].Draw()
                    icanvas = jZ + 2*nplots
                    canvas.cd(icanvas)
                    plotList[2].Draw("zcolpcol")
                    icanvas = jZ + 3*nplots
                    canvas.cd(icanvas)
                    plotList[3].Draw("zcolpcol")

            # set it so that python doesn't own these ROOT object
            for key in plotDict.keys():
                for plot in plotDict[key]:
                    SetOwnership(plot,False)

            # save canvas in the output Dictionary
            donutDict["canvas"] = canvas

        # all done
        return donutDict
            
            

    def calcHexapod(self,dictOfResults):

        zResultDict = dictOfResults["z4ResultDict"]
        z5ResultDict = dictOfResults["z5ResultDict"]
        z6ResultDict = dictOfResults["z6ResultDict"]
        z7ResultDict = dictOfResults["z7ResultDict"]
        z8ResultDict = dictOfResults["z8ResultDict"]

        # if any problems, print out error message
        try:

            # build aberration Column vector
            # average the redundant rotations of Astigmatism
            aveZern5ThetaX = 0.5 * (z5ResultDict["thetax"] + z6ResultDict["thetay"])
            aveZern6ThetaX = 0.5 * (z6ResultDict["thetax"] - z5ResultDict["thetay"])

            aveZern5ThetaXErr = numpy.sqrt( 0.25 * (z5ResultDict["thetaxErr"]*z5ResultDict["thetaxErr"] + z6ResultDict["thetayErr"]*z6ResultDict["thetayErr"]) )
            aveZern6ThetaXErr = numpy.sqrt( 0.25 * (z6ResultDict["thetaxErr"]*z6ResultDict["thetaxErr"] + z5ResultDict["thetayErr"]*z5ResultDict["thetayErr"]) )

            # build a weighted average
            try:
                wgtaveZern5ThetaXIErr2 =  ( 1.0/numpy.power(z5ResultDict["thetaxErr"],2) + 1.0/numpy.power(z6ResultDict["thetayErr"],2) )
                wgtaveZern5ThetaX = ( z5ResultDict["thetax"]/numpy.power(z5ResultDict["thetaxErr"],2) + z6ResultDict["thetay"]/numpy.power(z6ResultDict["thetayErr"],2) ) / wgtaveZern5ThetaXIErr2
                wgtaveZern5ThetaXErr = numpy.sqrt(1.0/wgtaveZern5ThetaXIErr2)
                wgtaveZern6ThetaXIErr2 =  ( 1.0/numpy.power(z6ResultDict["thetaxErr"],2) + 1.0/numpy.power(z5ResultDict["thetayErr"],2) )
                wgtaveZern6ThetaX = ( z6ResultDict["thetax"]/numpy.power(z6ResultDict["thetaxErr"],2) - z5ResultDict["thetay"]/numpy.power(z5ResultDict["thetayErr"],2) ) / wgtaveZern6ThetaXIErr2
                wgtaveZern6ThetaXErr = numpy.sqrt(1.0/wgtaveZern6ThetaXIErr2)
            except:
                wgtaveZern5ThetaX = aveZern5ThetaX
                wgtaveZern6ThetaX = aveZern6ThetaX
                wgtaveZern5ThetaXErr = aveZern5ThetaXErr
                wgtaveZern6ThetaXErr = aveZern6ThetaXErr
                
            # use the regular average here
            aberrationList = (aveZern5ThetaX,aveZern6ThetaX,z7ResultDict["delta"],z8ResultDict["delta"])
            aberrationColVec = numpy.matrix(aberrationList).transpose()
            
            # build aberration Error Column vector
            aberrationErrList = (aveZern5ThetaXErr,aveZern6ThetaXErr,z7ResultDict["deltaErr"],z8ResultDict["deltaErr"])
            aberrationErrArr = numpy.array(aberrationErrList)
            aberrationErrSquaredArr = aberrationErrArr * aberrationErrArr
            aberrationErrSquaredColVec = numpy.matrix(aberrationErrSquaredArr).transpose()
            
            # calculate hexapod vector        
            hexapodM = self.alignmentMatrix * aberrationColVec
            hexapodVector = hexapodM.A

            # calculated hexapod vector error
            hexapodErrSquaredM = self.alignmentMatrixSquared * aberrationErrSquaredColVec
            hexapodErrSquaredVector = hexapodErrSquaredM.A
            hexapodErrVector = numpy.sqrt(hexapodErrSquaredVector)
            
            # fill output dictionary (these are now in the Hexapod coordinate system)
            donutDict = {}
            
            donutDict["dodz"] = zResultDict["delta"] - self.paramDict["deltaZSetPoint"]
            donutDict["dodzErr"] = zResultDict["deltaErr"]
            
            donutDict["dodx"] = hexapodVector[0][0]
            donutDict["dody"] = hexapodVector[1][0]
            donutDict["doxt"] = hexapodVector[2][0]
            donutDict["doyt"] = hexapodVector[3][0]
        
            donutDict["dodxErr"] = hexapodErrVector[0][0]
            donutDict["dodyErr"] = hexapodErrVector[1][0]
            donutDict["doxtErr"] = hexapodErrVector[2][0]
            donutDict["doytErr"] = hexapodErrVector[3][0]

            # and for forward compatibility ALSO put this in a sub-dictionary called "donut_summary"
            # but with the MINUS SIGN to make these corrections instead of measurements
            # fill all DB variables into donut_summary also AJR 12/1/2012
            donutSummary = {}

            donutSummary["dodz"] = -(zResultDict["delta"] - self.paramDict["deltaZSetPoint"])
            donutSummary["dodzerr"] = zResultDict["deltaErr"]
            
            donutSummary["dodx"] = -donutDict["dodx"]
            donutSummary["dody"] = -donutDict["dody"]
            donutSummary["doxt"] = -donutDict["doxt"]
            donutSummary["doyt"] = -donutDict["doyt"]
        
            donutSummary["dodxerr"] = donutDict["dodxErr"]
            donutSummary["dodyerr"] = donutDict["dodyErr"]
            donutSummary["doxterr"] = donutDict["doxtErr"]
            donutSummary["doyterr"] = donutDict["doytErr"]
            
            # also fill all result variables into donutSummary
            # only donutSummary values are used in the downstream analysis, csv file and TTree
            # --- changing zdelta to z4delta (not backwards compatible!)

            infoList = ["delta","thetax","thetay","deltaErr","thetaxErr","thetayErr","meanDeltaBefore","rmsDeltaBefore","meanDeltaAfter","rmsDeltaAfter"]
            for key in dictOfResults:
                if key[0] == "z":
                    for datum in infoList:
                        datumName = "%s%s" % (key.replace("ResultDict",""),datum)
                        donutSummary[datumName] = dictOfResults[key][datum]

            # get number of donuts ultimately used
            ndonuts_used = len(zResultDict["deltaArrayAfter"])
            donutSummary["ndonuts_used"] = ndonuts_used

            # store in the output dictionary
            donutDict["donut_summary"] = donutSummary

        except:
            print "donutana: Not enough information for calcHexapod"
            donutDict = {}

        return donutDict
        
        
        

    def fitToRefMesh(self,refMesh,otherMesh,angleconversion=1.0):
        """ fit a mesh to a Reference mesh, return the fit parameters
        and errors and convenient arrays for plots
        """

        # resultsDict contents
        #      delta fit result
        #      thetax fit result
        #      thetay fit result
        #      error on delta fit result
        #      error on thetax fit result
        #      error on thetay fit result
        #      mean of differences before fit
        #      sigma of residual distribution
        #      array of differences before the fit
        #      array of differences after the fit

        resultDict = {}
                                     
        # fit the two meshes
        # note that this fits to 
        #       zDiff = otherMesh - interp of refMesh
        results,zDiff,xColumn,yColumn = refMesh.fitMesh(otherMesh,"RLM")

        if results!=None:
            # fill resultsDict
            resultDict["delta"] = results.params[2]
            resultDict["deltaErr"] = results.bse[2]
            
            # convert from [Values/Dictance] to another unit
            resultDict["thetax"] = results.params[0]*angleconversion
            resultDict["thetay"] = results.params[1]*angleconversion
            resultDict["thetaxErr"] = results.bse[0]*angleconversion
            resultDict["thetayErr"] = results.bse[1]*angleconversion

            # mean delta before the fit - no truncation now
            resultDict["meanDeltaBefore"] = scipy.mean(zDiff)
            resultDict["rmsDeltaBefore"] = scipy.std(zDiff)
            
            # rms of Differences after fit
            resultDict["meanDeltaAfter"] = scipy.mean(results.resid)
            resultDict["rmsDeltaAfter"] = scipy.std(results.resid)

            # array of differences before the fit
            resultDict["deltaArrayBefore"] = zDiff            

            # array of differences after the fit
            resultDict["deltaArrayAfter"] = results.resid
            resultDict["deltaArrayX"] = xColumn
            resultDict["deltaArrayY"] = yColumn
            
        return resultDict


    def calcMeshWgts(self,dictOfMeshes,type="nMAD",kNN=200,nMADCut=4.0):
        """ utility routine to recalculate weights for each mesh entry 
        calculation is done in-place
        """
        
        outDict = {}
        # calculate nMAD over kNN for each point in each mesh
        # then calculate a new weight: =1 for nMAD<cut and =0 for nMAD>cut
        # try removing tails 
        if type=="nMAD":
            plotList = []
            for iZ in range(4,15+1):
                key = "z%dMesh" % (iZ)
                if dictOfMeshes.has_key(key):
                    thisMesh = dictOfMeshes[key]
                    nMADArr = thisMesh.calcWgtnMAD(kNN=kNN,nMADCut=nMADCut)

                    hnMAD = hfillhist("z%dnMAD" % (iZ),"z%d nMAD Distribution" % (iZ),nMADArr,100,-10.,10.)
                    plotList.append(hnMAD)

            # unique name for our canvas
            tstr = "canvas" + str(time.time())
            nplots = len(plotList)
            nCols = 4
            nRows = (nplots-1)/nCols + 1
            canvas = TCanvas(tstr,tstr,300*nRows,300*nCols)
            canvas.Divide(nRows,nCols)
            
            for i in range(nplots):
                canvas.cd(i+1)
                plotList[i].Draw()
                SetOwnership(plotList[i],False)

            outDict["canvas"] = canvas

        else:
            print "donutana Error: bad type ",type," in calcMeshWgts"

        return outDict


    def cullAllMeshes(self,dictOfMeshes,cullCut=0.01):
        """ cull a dictionary of Meshes, where all Meshes are of the same donuts 
        calculates the product of weights for each donut, and removes donuts from all meshes which fail a 
        cut on the weight
        """

        # loop over the Coords and coalesce the Wgts of the points

        aMesh = dictOfMeshes[dictOfMeshes.keys()[0]]   # this is just the first mesh
        dictOfWgts = {}
        for coord in aMesh.coordList:

            if aMesh.pointsArray.has_key(coord):
                npts = aMesh.pointsArray[coord].shape[0]
                if npts>0 :
                    thisWgt = numpy.ones(npts)
                    # combine weights with the product 
                    for mesh in dictOfMeshes:
                        thisMesh = dictOfMeshes[mesh]
                        thisWgt = thisWgt * thisMesh.pointsArray[coord][:,3]

                    # now insert this weight for all meshes
                    for mesh in dictOfMeshes:
                        thisMesh = dictOfMeshes[mesh]
                        thisMesh.pointsArray[coord][:,3] = thisWgt

        # finally cull each mesh!
        for mesh in dictOfMeshes:
            thisMesh = dictOfMeshes[mesh]
            thisMesh.cullMesh(cullCut)

        return 0



    def analyzeRzero(self,rzeroMesh,chi2Mesh,neleMesh):
        """ utility routine to analyze fitted rzero and/or chi2 vs. nele
        """
        resultDict = {}

        # only works for FandAOnly
        if self.paramDict["sensorSet"] == "FandAOnly" :

            extraCoords = ["FS1","FS3","FN2","FN4"]
            intraCoords = ["FS2","FS4","FN1","FN3"]
            bothCoords = {"extra":extraCoords,"intra":intraCoords}
            # see Evernote "Seeing from r0 fits or Chi2"
            const = {"extra":-2.891,"intra":-2.850}
            slope = {"extra":1.215,"intra":1.191}

            for focal in bothCoords.keys():
                coords = bothCoords[focal]
                rzeroAllPntsL = []
                neleAllPntsL = []
                chi2AllPntsL = []
                for aCoord in coords:
                    rzeroArr = rzeroMesh.pointsArray[aCoord]
                    if len(rzeroArr)>0:
                        rzeroAllPntsL.extend(rzeroArr[:,2].tolist())

                    chi2Arr = chi2Mesh.pointsArray[aCoord]
                    if len(chi2Arr)>0:
                        chi2AllPntsL.extend(chi2Arr[:,2].tolist())

                    neleArr = neleMesh.pointsArray[aCoord]
                    if len(neleArr)>0:
                        neleAllPntsL.extend(neleArr[:,2].tolist())

                rzeroPnts = numpy.array(rzeroAllPntsL)
                chi2Pnts = numpy.array(chi2AllPntsL)
                nelePnts = numpy.array(neleAllPntsL)

                # calculate median and MAD of rzero
                name = "rzeroMedian" + focal
                rzeroMedian = numpy.median(rzeroPnts)
                resultDict[name] = rzeroMedian
                name = "rzeroMAD" + focal
                resultDict[name] = numpy.median(numpy.abs(rzeroPnts-rzeroMedian))

                # calculate median and MAD of corrected Chi2
                corrchiPnts = -(numpy.log10(chi2Pnts) - (const[focal] + slope[focal]*numpy.log10(nelePnts)))
                name = "corrchiMedian" + focal
                corrchiMedian = numpy.median(corrchiPnts)
                resultDict[name] = corrchiMedian
                name = "corrchiMAD" + focal
                resultDict[name] = numpy.median(numpy.abs(corrchiPnts-corrchiMedian))

        return resultDict
 
