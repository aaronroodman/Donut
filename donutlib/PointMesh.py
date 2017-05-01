#
# Python Class to encapsulate a 2-D mesh of points
# multiple "Coordinate" systems are possible, each with separate
# interpolation grids
#
import numpy
import scipy
import pandas as pd
from scipy.interpolate import SmoothBivariateSpline
from scipy.interpolate.rbf import Rbf
from scipy.interpolate import griddata
from scipy.spatial import cKDTree
from donutlib.IDWInterp import IDWInterp
from donutlib.decamutil import decaminfo
import bisect

try:
    from scipy import stats
except:
    print("PointMesh: Could not import scipy.stats")
    
import numpy.lib.index_tricks as itricks
import statsmodels.api as sm
try:
    from matplotlib import cm,colors
    from mpl_toolkits.mplot3d import axes3d
    from matplotlib import pyplot as plt
except:
    print("PointMesh: Could not load matplotlib")
    
try:
    from ROOT import TStyle,gStyle,TCanvas,TGraph2D
except:
    print("PointMesh: Could not load ROOT")

import pdb


class PointMesh(object):
    """ PointMesh encapsulates a 2-D mesh of points.  The points have y,x,index
    dimensions.  Each value of index may provide its own coordinate system, or there may be a
    single coordinate system. Multiple interpolation methods are implemented.

    Aaron Roodman (C) SLAC National Accelerator Laboratory, Stanford University 2012.
    """
    def readPointsFromFile(self,fileName):

        # read in pointsArray from file, use pandas now
        dataPoints = pd.read_csv(fileName, delim_whitespace=True,
                    header=None,
                    dtype={'Sensor': '|S3', 'x': numpy.float64, 'y': numpy.float64,
                           'z': numpy.float64, 'w': numpy.float64},
                    names=['Sensor', 'x', 'y', 'z', 'w'])

        # convert the raw pandas format to a dictionary keyed by sensor name
        # containing a numpy array of x,y,z,w
        self.pointsArray = {}
        dataPoints['Sensor'] = dataPoints['Sensor'].str.decode('utf-8')
        for iCoord in self.coordList:
            selone = (dataPoints['Sensor']==iCoord)
            xarr = numpy.array(dataPoints[selone]['x'].tolist())
            yarr = numpy.array(dataPoints[selone]['y'].tolist())
            zarr = numpy.array(dataPoints[selone]['z'].tolist())
            warr = numpy.array(dataPoints[selone]['w'].tolist())
            arr = numpy.dstack([xarr,yarr,zarr,warr])[0]   #also numpy.array(zip( ))  works too
            self.pointsArray[iCoord] = arr

    def checkMethod(self,myMethod,methodVal=None):
        if self.myMethod=='sbs' or self.myMethod=='rbf' or self.myMethod=='tmean' or self.myMethod=='grid' or self.myMethod == "idw" or self.myMethod == "bmedian":
            ok = True
        else:
            ok = False
        return ok            

    def __init__(self,coordList,gridArray,pointsArray=None,pointsFile=None,myMethod='sbs',methodVal=None,debugFlag=False,title=""):
        """ initialize the class
        """
        self.debugFlag = debugFlag

        self.title = title

        self.nCoord = len(coordList)
        self.coordList = coordList

        # vignetting  HARDCODED
        self.radiusVignetted = 225.0

        # decam info
        self.decam = decaminfo()

        # Array with interpolation grid size and number of points
        # should contain a 4 tuple of [ny,ylo,yhi,nx,xlo,xhi] for each nCoord
        # these ylo,yhi and xlo,xhi should be the region EDGES
        self.gridArray = gridArray.copy()
        
        # contains an array of Points for each nCoord
        # [npoints,0] is the X value
        # [npoints,1] is the Y value
        # [npoints,2] is the Z value
        # [npoints,3] is the Weight value
        # the arrays are stored in a Python dictionary, keyed by nCoord index
        # NOTE:  nothing in this code enforces that pointArray has this content!!!
        if pointsArray!=None:
            self.pointsArray = pointsArray.copy()
        elif pointsFile!=None:
            self.readPointsFromFile(pointsFile)
        else:
            print("PointMesh: no input points")

        # construct the interpolation grid for each coordinate system
        self.myMethod = myMethod
        self.methodVal = methodVal
        self.interpPresent = False
        if self.checkMethod(myMethod,methodVal):
            self.interpPresent = True
            self.makeInterpolation(self.myMethod,self.methodVal)

    def __getstate__(self):
        stateDict = {}
        keysToPickle = ['nCoord','gridArray','radiusVignetted','pointsArray','myMethod','coordList','debugFlag','methodVal','title','decam']
        for key in keysToPickle:
            stateDict[key] = self.__dict__[key]
        return stateDict

    def __setstate__(self,state):
        for key in state:
            self.__dict__[key] = state[key]

        self.interpPresent = False                
        if self.checkMethod(self.myMethod,self.methodVal):
            self.interpPresent = True
            self.makeInterpolation(self.myMethod,self.methodVal)
        


    def makeGrid(self,ny,ylo,yhi,nx,xlo,xhi):
        # make a grid using mgrid with grid edges from xlo->xhi and ylo->yhi

        xdelta = (xhi-xlo)/nx
        xlolimit = xlo + 0.5*xdelta
        xhilimit = xhi - 0.5*xdelta
        xusen = nx 

        ydelta = (yhi-ylo)/ny
        ylolimit = ylo + 0.5*ydelta
        yhilimit = yhi - 0.5*ydelta
        yusen = ny 

        # these edges return a grid:
        # [-3840., -3328., -2816., -2304., -1792., -1280.,  -768.,  -256.,
        #    256.,   768.,  1280.,  1792.,  2304.,  2816.,  3328.,  3840.]
        # when used with lo=-4096. and hi=4096. and n=16
        yGrid,xGrid = itricks.mgrid[ylolimit:yhilimit:1j*yusen,xlolimit:xhilimit:1j*xusen]

        # now ALSO make a grid of size ny+1,nx+1 where the values are the grid edges
        yEdge,xEdge = itricks.mgrid[ylo:yhi:1j*(yusen+1),xlo:xhi:1j*(xusen+1)]

        return yGrid,xGrid,yEdge,xEdge

    def makeInterpolation(self,myMethod="sbs",methodVal=None):
        # construct the interpolation grids for all coordinate systems

        # dictionaries with interpolation grids and coordinate systems
        self.interpGrids = {}
        self.interpEdges = {}
        self.interpValues = {}
        self.interpSpline = {}
        self.interpRbf = {}
        self.interpIDW = {}
        self.interpBMedian = {}
        self.interpTMean = {}
        self.interpTStd = {}

        # loop over Coordinate systems
        for iCoord in self.coordList:

            # build cell-centers for the interpolation grid
            ny,ylo,yhi,nx,xlo,xhi = self.gridArray[iCoord]
            yGrid,xGrid,yEdge,xEdge = self.makeGrid(ny,ylo,yhi,nx,xlo,xhi)
            self.interpGrids[iCoord] = [xGrid,yGrid]
            self.interpEdges[iCoord] = [xEdge,yEdge]

            data = self.pointsArray[iCoord]
            if self.debugFlag:
                print("PointMesh: At ",iCoord,"we have ",data.shape[0]," points")
            npts = data.shape[0]

            # check number of points
            if npts>=5:
                xData = data[:,0] 
                yData = data[:,1] 
                zData = data[:,2]

                if myMethod =="sbs":

                    # SmoothBivariateSpline
                    if npts>600:
                        self.interpSpline[iCoord] = SmoothBivariateSpline(xData,yData,zData,bbox=[xlo,xhi,ylo,yhi],kx=4,ky=4,s=1.e6)
                    elif npts>=100:
                        self.interpSpline[iCoord] = SmoothBivariateSpline(xData,yData,zData,bbox=[xlo,xhi,ylo,yhi],kx=3,ky=3,s=1.e6)
                    elif npts>9:
                        self.interpSpline[iCoord] = SmoothBivariateSpline(xData,yData,zData,bbox=[xlo,xhi,ylo,yhi],kx=2,ky=2,s=1.e6)
                    else:
                        self.interpSpline[iCoord] = SmoothBivariateSpline(xData,yData,zData,bbox=[xlo,xhi,ylo,yhi],kx=1,ky=1,s=1.e7)
                    self.interpValues[iCoord] = self.interpSpline[iCoord].ev(xGrid.reshape((ny*nx)),yGrid.reshape((ny*nx))).reshape((ny,nx))


                elif myMethod =="rbf":

                    self.interpRbf[iCoord] = Rbf(xData,yData,zData)
                    self.interpValues[iCoord] = self.interpRbf[iCoord](xGrid.reshape((ny*nx)),yGrid.reshape((ny*nx))).reshape((ny,nx))

                elif myMethod =="tmean":

                    # use the truncated mean for each Coord -- very very simple!!
                    zstd = stats.tstd(zData)
                    zmean = stats.tmean(zData)
                    ztmean = stats.tmean(zData,(zmean-3.*zstd,zmean+3.*zstd))
                    ztstd = stats.tstd(zData,(zmean-3.*zstd,zmean+3.*zstd))
                    self.interpTMean[iCoord] = ztmean
                    self.interpTStd[iCoord] = ztstd
                    self.interpValues[iCoord] = ztmean*numpy.ones((ny,nx))

                elif myMethod =="bmedian":

                    # use the median for each bin in each Coord
                    self.interpBMedian[iCoord] = numpy.zeros((ny,nx))
                    zvalL = []
                    xbin = numpy.digitize(xData,xEdge[0,:]) - 1
                    ybin = numpy.digitize(yData,yEdge[:,0]) - 1
                    for i in range(nx):
                        for j in range(ny):
                            ok = numpy.logical_and.reduce((xbin==i,ybin==j))
                            zHere = zData[ok]
                            nEntry = zHere.shape[0]
                            if nEntry>=1:
                                self.interpBMedian[iCoord][j,i] = numpy.median(zHere)
                            else:
                                # need to do something better!
                                self.interpBMedian[iCoord][j,i] =0.
                    # fill interpValues, need to match order of locations in xGrid and yGrid
                    self.interpValues[iCoord] = self.interpBMedian[iCoord].copy()          

                elif myMethod == "grid":
                    Z = griddata((xData,yData),zData,(xGrid.reshape((ny*nx)),yGrid.reshape((ny*nx))),"linear")
                    self.interpValues[iCoord] =  Z.reshape(xGrid.shape)

                elif myMethod == "idw":
                    # July 15, 2013 - change to use epsilon=1.0 (mm) to set a cutoff in the distance
                    # this will make small changes in the results for all Donuts
                    if methodVal != None :
                        usekNN = methodVal[0]
                        useEpsilon = methodVal[1]
                    else:
                        usekNN = 4
                        useEpsilon = 1.0
                    self.interpIDW[iCoord] = IDWInterp(xData,yData,zData,kNN=usekNN,epsilon=useEpsilon)
                    self.interpValues[iCoord] = self.interpIDW[iCoord].ev(xGrid.reshape((ny*nx)),yGrid.reshape((ny*nx))).reshape((ny,nx))
                    
                else:
                    self.interpSpline[iCoord] = None
                    self.interpValues[iCoord] = numpy.zeros(xGrid.shape)

            else:
                self.interpTMean[iCoord] = 0.0
                self.interpBMedian[iCoord] = None
                self.interpTStd[iCoord] = 0.0
                self.interpSpline[iCoord] = None
                self.interpRbf[iCoord] = None
                self.interpValues[iCoord] = numpy.zeros(xGrid.shape)



    def doInterp(self,iCoord,x,y):
        """ get interpolated value for iCoord,x,y
        """

        # put all the returns at the bottom - to check that Z is not a "nan"
        
        if type(iCoord) == int:
            iCoord = str(iCoord)
            
        if self.myMethod =="sbs":
            if self.interpSpline[iCoord] != None:
                return self.interpSpline[iCoord].ev(x,y)
            else:
                Z = numpy.zeros(x.shape)
                return Z
            
        elif self.myMethod == "rbf":
            if self.interpRbf[iCoord] != None:
                return self.interpRbf[iCoord](x,y)
            else:
                Z = numpy.zeros(x.shape)
                return Z

        elif self.myMethod == "idw":
            if iCoord in self.interpIDW and self.interpIDW[iCoord] != None:
                return self.interpIDW[iCoord].ev(x,y)
            else:
                Z = numpy.zeros(x.shape)
                return Z

        elif self.myMethod == "tmean":
            if self.interpTMean[iCoord] != None:
                return self.interpTMean[iCoord]*numpy.ones(x.shape)
            else:
                Z = numpy.zeros(x.shape)
                return Z

        elif self.myMethod == "bmedian":
            if self.interpBMedian[iCoord] != None:
                xEdge,yEdge = self.interpEdges[iCoord]
                xbin = numpy.digitize(x,xEdge[0,:]) - 1
                ybin = numpy.digitize(y,yEdge[:,0]) - 1
                return self.interpBMedian[iCoord][ybin,xbin]
            else:
                Z = numpy.zeros(x.shape)
                return Z

        elif self.myMethod == "grid":
            # unfortunately griddata needs to be filled with the points each time
            Z = griddata((xData,yData),zData,(x,y),"linear")
            return Z
        
        else:
            Z = numpy.zeros(x.shape)
            return Z

            
        

    def fitMesh(self,otherMesh,method="OLS"):
        """ fit another Mesh to this one, find offset,tip,tilt to match the 2 meshes
        compare self's Interpolation mesh and the discrete points from otherMesh
        note that we fit to zDiff = otherMesh - self.interp
        """

        interpData = {}
        npoints = 0
        # loop over coordinate systems
        for iCoord in self.coordList:

            # loop over points in otherMesh
            if iCoord in otherMesh.pointsArray:
                otherData = otherMesh.pointsArray[iCoord]
                if otherData.shape[0]>0:
                    xOtherData = otherData[:,0]
                    yOtherData = otherData[:,1]
                    npoints = npoints + otherData.shape[0]
                    # interpData[iCoord] = self.interpSpline[iCoord].ev(xOtherData,yOtherData)
                    interpData[iCoord] = self.doInterp(iCoord,xOtherData,yOtherData)

        # now we can compare zOtherData and interpData

        # difference vector
        zDiff = numpy.zeros(npoints)

        # columns with 1,x,y
        oneColumn = numpy.ones(npoints)
        xColumn = numpy.zeros(npoints)
        yColumn = numpy.zeros(npoints)

        nstart = 0
        for iCoord in self.coordList:
            if iCoord in otherMesh.pointsArray:
                npts = otherMesh.pointsArray[iCoord].shape[0]
                if npts>0 :
                    zDiff[nstart:nstart+npts] = otherMesh.pointsArray[iCoord][:,2] - interpData[iCoord]
                    xColumn[nstart:nstart+npts] = otherMesh.pointsArray[iCoord][:,0]
                    yColumn[nstart:nstart+npts] = otherMesh.pointsArray[iCoord][:,1]
                    nstart = nstart + npts

        # build matrix
        # formula to fit is ThetaX*y + ThetaY*x + DeltaZ
        # so note that ThetaX,ThetaY are in units of [Value/Distance] and DeltaZ are in units of [Value]
        aMatrix = numpy.vstack([yColumn,xColumn,oneColumn]).T

        # solve linear least square
        #thetax,thetay,deltaz = numpy.linalg.lstsq(aMatrix, zDiff)[0]

        # make sure there is something to fit
        if nstart>3:
            # use Regression method or use sm.RLM for robust fitting
            if method=="OLS":
                linearModel = sm.OLS(zDiff,aMatrix)
            elif method=="RLM":
                linearModel = sm.RLM(zDiff,aMatrix)
            
            try:
                results = linearModel.fit()
            except:
                print("PointMesh.py:  fit failed")
                results = None
                
            if self.debugFlag:
                print("PrintMesh: ThetaX = ",results.params[0]," +- ",results.bse[0])
                print("PrintMesh: ThetaY = ",results.params[1]," +- ",results.bse[1])
                print("PrintMesh: DeltaZ = ",results.params[2]," +- ",results.bse[2])

            if method=="RLM" and results!=None :
                mean = (results.resid*results.weights).sum()/results.weights.sum()
                rms = numpy.sqrt( (numpy.power(results.resid-mean,2)*results.weights).sum() / results.weights.sum()  )
                if self.debugFlag:
                    print("PointMesh: Fraction of points used = ",results.weights.sum()/results.nobs)
                    print("PointMesh: Mean, RMS with weights = ",mean,rms)

            # load results.weights into otherMesh's weights
            if results!=None:
                nstart = 0
                for iCoord in otherMesh.coordList:
                    if iCoord in otherMesh.pointsArray:
                        npts = otherMesh.pointsArray[iCoord].shape[0]
                        if npts>0:
                            otherMesh.pointsArray[iCoord][:,3] = results.weights[nstart:nstart+npts]
                            nstart = nstart + npts
                    
        else:
            results = None

        return results,zDiff,xColumn,yColumn

    def diffMesh(self,otherMesh):
        # find average of another Mesh to this one, 

        interpValues = {}
        npoints = 0
        # loop over coordinate systems
        for iCoord in self.coordList:

            # loop over points in otherMesh
            otherData = otherMesh.pointsArray[iCoord]
            if otherData.shape[0]>0:
                xOtherData = otherData[:,0]
                yOtherData = otherData[:,1]
                npoints = npoints + otherData.shape[0]
#                interpValues[iCoord] = self.interpSpline[iCoord].ev(xOtherData,yOtherData)
                interpValues[iCoord] = self.doInterp(iCoord,xOtherData,yOtherData)

        # now we can compare zOtherData and interpValues

        # difference vector
        zDiff = numpy.zeros(npoints)

        # columns with 1,x,y
        oneColumn = numpy.ones(npoints)
        xColumn = numpy.zeros(npoints)
        yColumn = numpy.zeros(npoints)

        nstart = 0
        for iCoord in self.coordList:
            npts = otherMesh.pointsArray[iCoord].shape[0]
            if npts>0:
                zDiff[nstart:nstart+npts] = otherMesh.pointsArray[iCoord][:,2] - interpValues[iCoord]
                xColumn[nstart:nstart+npts] = otherMesh.pointsArray[iCoord][:,0]
                yColumn[nstart:nstart+npts] = otherMesh.pointsArray[iCoord][:,1]
                nstart = nstart + npts
                
        return zDiff,xColumn,yColumn




#    def diffMeshInterp(self,otherMesh):
#        """ find difference of another Mesh to this one, using interpolation for both
#        """
#
#        # loop over coordinate systems
#        for iCoord in self.coordList:
#
#            xGrid,yGrid = self.interpGrids[iCoord]
#
#            # loop over this meshes interpolation Grid
#            # then interpolate the other other Mesh on this Grid too
#
#                interpValues[iCoord] = self.doInterp(iCoord,xOtherData,yOtherData)
#
#        # return diff,x,y, arrays
#        return zDiff,xColumn,yColumn








    def singleToMultiple(self,input):
        # utility to convert a single flat array to multiple arrays in a Dictionary, one for each "Coordinate"
        # given an input with the same dimensions as the pointsArray
        nstart = 0
        multArray = {}
        for iCoord in self.coordList:
            aArray = self.pointsArray[iCoord]
            npoints = aArray.shape[0]
            if npoints>0:
                tempArray = input[nstart:nstart+npoints]
                nstart = nstart + npoints
                multArray[iCoord] = tempArray
                
        return multArray
 

    def getXYZpoints(self,coordList=None):
        """ utility to return separate X,Y,Z arrays for all input points
        """

        if coordList==None:
            coordList = self.coordList

        dataX = []
        dataY = []
        dataZ = []
        
        # loop over Coordinate systems
        for iCoord in coordList:

            data = self.pointsArray[iCoord]
            npts = data.shape[0]
            if npts>0:
                dataX.extend(data[:,0].tolist())
                dataY.extend(data[:,1].tolist())
                dataZ.extend(data[:,2].tolist())

        return numpy.array(dataX),numpy.array(dataY),numpy.array(dataZ)


    def getXYZWpoints(self,coordList=None):
        """ utility to return separate X,Y,Z,W arrays for all input points
        """
        if coordList==None:
            coordList = self.coordList
        
        dataX = []
        dataY = []
        dataZ = []
        dataW = []
        
        # loop over Coordinate systems
        for iCoord in coordList:

            data = self.pointsArray[iCoord]
            npts = data.shape[0]
            if npts>0:
                dataX.extend(data[:,0].tolist())
                dataY.extend(data[:,1].tolist())
                dataZ.extend(data[:,2].tolist())
                dataW.extend(data[:,3].tolist())

        return numpy.array(dataX),numpy.array(dataY),numpy.array(dataZ),numpy.array(dataW)

    def autoRange(self,zmin=None,zmax=None,coordList=None):
        """ determine X,Y,Z limits
        """

        if coordList==None:
            coordList = self.coordList

        # autorange in X,Y
        minX = 1.e50
        maxX = -1.e50
        minY = 1.e50
        maxY = -1.e50
        minZ = 1.e50
        maxZ = -1.e50
        for iCoord in coordList:

            #if self.interpPresent:
            #    X,Y = self.interpGrids[iCoord]
            #    Z = self.interpValues[iCoord]
            #else:
            X,Y,Z = self.getXYZpoints(coordList=coordList)

            minX = numpy.minimum(minX,X.min())
            maxX = numpy.maximum(maxX,X.max())
            minY = numpy.minimum(minY,Y.min())
            maxY = numpy.maximum(maxY,Y.max())
            minZ = numpy.minimum(minZ,Z.min())
            maxZ = numpy.maximum(maxZ,Z.max())

        if zmin!=None:
            minZ = zmin
        if zmax!=None:
            maxZ = zmax

        return minX,maxX,minY,maxY,minZ,maxZ


    def plotDataROOT(self,name='plotDataROOT',ztitle='Z',zmin=None,zmax=None):

        self.canvasData = TCanvas("canvasData","canvasData",0,0,300,300)
        gStyle.SetOptStat(0)

        minX,maxX,minY,maxY,minZ,maxZ = self.autoRange(zmin,zmax)

        X,Y,Z = self.getXYZpoints()
        n = X.shape[0]

        self.graphData = TGraph2D(name,ztitle,n,X,Y,Z)

        self.graphData.SetMaximum(maxZ)
        self.graphData.SetMinimum(minZ)
        gStyle.SetHistTopMargin(0)
        self.graphData.Draw("zpcol")
        return self.graphData


    def plotdMi(self,title='Z',range=None,coordList=None,bins=100):
        """ plot the difference between the data and the interpolation
        """
        # get coordList
        if coordList == None:
            coordList = self.coordList
        
        # reserve a Figure and 3D axes
        fig = plt.figure(figsize=(10,7.5))

        # set ranges before too
        self.ax.set_xlabel('Data - Interpolation')
        self.ax.set_xlim(minX, maxX)

        X,Y,Z = self.getXYZpoints(coordList)
        Zi = self.doInterp(X,Y) 
        plt.hist(Z-Zi,nbins=nbins,range=range,title=title)

        # plot it!
        plt.show()

        return fig

    

    # to fix current color bug
    def forceUpdate(self,event):
        self.colorScatterPlot.changed()
            

    def plotMeshMPL(self,same=False,plotInterp=True,plotColor=True,plotData=True,ztitle='Z',zmin=None,zmax=None,nstride=1,title="",coordList=None):
        """ plot the Mesh
        """
        # get coordList
        if coordList == None:
            coordList = self.coordList
        
        # reserve a Figure and 3D axes
        if not same:
            self.fig = plt.figure(figsize=(10,7.5))
            self.ax = self.fig.gca(projection='3d')

        # autorange in X,Y
        minX,maxX,minY,maxY,minZ,maxZ = self.autoRange(zmin,zmax,coordList=coordList)

        # set ranges before too
        self.ax.set_xlabel('X')
        self.ax.set_xlim(minX, maxX)
        self.ax.set_ylabel('Y')
        self.ax.set_ylim(minY, maxY)
        self.ax.set_zlabel(ztitle)
        self.ax.set_zlim(minZ,maxZ)    
        anorm = colors.Normalize(minZ,maxZ,True)

        if plotData:
            X,Y,Z = self.getXYZpoints(coordList)
            self.colorScatterPlot = self.ax.scatter(X, Y, Z, c=Z, s=25., marker='o',edgecolors='none',cmap=cm.jet,norm=anorm)
            self.fig.colorbar(self.colorScatterPlot)
            self.fig.canvas.mpl_connect('draw_event',self.forceUpdate)


        # plot Interpolation grids
        if self.interpPresent and (plotColor or plotInterp):

            # loop over Extensions
            for iCoord in coordList:
                
                X,Y = self.interpGrids[iCoord]
                R = numpy.sqrt(X*X + Y*Y)
                XVig = numpy.ma.masked_where(R>self.radiusVignetted,X)
                YVig = numpy.ma.masked_where(R>self.radiusVignetted,Y)
                Z = self.interpValues[iCoord]
                ZVig = numpy.ma.masked_where(R>self.radiusVignetted,Z)

                if plotColor:

                    cset = self.ax.contourf(XVig, YVig, ZVig, zdir='z', offset=minZ, cmap=cm.jet, norm=anorm)
                if plotInterp:

                    #ax.plot_surface(X,Y,Z, rstride=nstride, cstride=nstride, alpha=0.3)
                    #ax.plot_surface(X,Y,Z, rstride=nstride, cstride=nstride, cmap=cm.jet)
                    self.ax.plot_wireframe(XVig,YVig,ZVig, rstride=nstride, cstride=nstride, color='0.25')


        # show axes the way I like
        self.ax.view_init(azim=-120)

        # set ranges last
        self.ax.set_xlabel('X')
        self.ax.set_xlim(minX, maxX)
        self.ax.set_ylabel('Y')
        self.ax.set_ylim(minY, maxY)
        self.ax.set_zlabel(ztitle)
        self.ax.set_zlim(minZ,maxZ)
        self.ax.set_title(title)

        # plot it!
        plt.show()

        return self.fig,self.ax


    def format_coord(self, x, y):
        """ given x,y find z value from the interpolation grid
        """
        extname = self.decam.getSensor(x,y)
        if extname==None:
            return     'x=%1.4f, y=%1.4f               '%(x, y)
        else:
            if extname in self.interpEdges:
                xEdges = self.interpEdges[extname][0][0,:]
                yEdges = self.interpEdges[extname][1][:,0]
            else:
                return     'x=%1.4f, y=%1.4f               '%(x, y)

            # get the right index into the interpValues array
            indY = bisect.bisect_left(yEdges,y)-1
            indX = bisect.bisect_left(xEdges,x)-1

        if extname in self.interpValues:
            shape = self.interpValues[extname].shape
            if indX>=0 and indX<shape[1] and indY>=0 and indY<shape[0]:
                z = self.interpValues[extname][indY,indX]
                return 'x=%1.4f, y=%1.4f, z=%1.4f, %s'%(x, y, z,extname)
            else:
                print("problem: ",extname,x,y)
                return 'x=%1.4f, y=%1.4f               '%(x, y)
        else:
            return     'x=%1.4f, y=%1.4f               '%(x, y)

            



    def plotMeshMPL2D(self,zmin=None,zmax=None,title="",coordList=None,overlay=False,overlayC=None,cmap=cm.jet):
        """ plot the Mesh
        """

        # get coordList
        if coordList == None:
            coordList = self.coordList
        
        # reserve a Figure and 3D axes
        if overlay:
            # use the current figure and axes
            self.fig = plt.gcf()
            self.ax = plt.gca()
        else:
            self.fig = plt.figure(figsize=(10,7.5))
            self.ax = self.fig.gca()

        # autorange in X,Y
        if overlay:
            minZ,maxZ = overlayC.get_clim()
            minX,maxX = self.ax.get_xlim()
            minY,maxY = self.ax.get_ylim()
        else:
            minX,maxX,minY,maxY,minZ,maxZ = self.autoRange(zmin,zmax,coordList=coordList)

        # hard code X,Y limits to be DECam boundaries
        minX = -235.
        maxX = 235.
        minY = -235.
        maxY = 235.
        
        # set ranges before too
        if not overlay:
            self.ax.set_xlabel('X',fontsize=24)
            self.ax.set_xlim(minX, maxX)
            self.ax.set_ylabel('Y',fontsize=24)
            self.ax.set_ylim(minY, maxY)
            
        anorm = colors.Normalize(minZ,maxZ,True)

        # loop over Extensions
        for iCoord in coordList:
                
            # pcolor needs X,Y values of the edges of the bins...
            X,Y = self.interpEdges[iCoord]
            Xcen,Ycen = self.interpGrids[iCoord]
            Rcen = numpy.sqrt(Xcen*Xcen + Ycen*Ycen)
            #R = numpy.sqrt(X*X + Y*Y)
            #XVig = numpy.ma.masked_where(R>self.radiusVignetted,X)
            #YVig = numpy.ma.masked_where(R>self.radiusVignetted,Y)
            Z = self.interpValues[iCoord]
            ZVig = numpy.ma.masked_where(Rcen>self.radiusVignetted,Z)

            self.cset = self.ax.pcolor(X, Y, ZVig, cmap=cmap, norm=anorm)
                
        # set ranges last
        if not overlay:
            self.ax.set_xlabel('X',fontsize=24)
            self.ax.set_xlim(minX, maxX)
            self.ax.set_ylabel('Y',fontsize=24)
            self.ax.set_ylim(minY, maxY)
            if title=="":
                title = self.title
            self.ax.set_title(title,fontsize=24)

            # add the colorbar
            self.fig.colorbar(self.cset)

            # setup plot so that interactively you can use the mouse to show the Z value!        
            self.ax.format_coord = self.format_coord
                    

        # plot it!
        plt.show()

        return self.fig,self.ax,self.cset



    def calcWgtDtoInterp(self,maxdist=0.1):
        # set the weight for the mesh by the distance from the actual value to the interpolated value
        # assuming that the interpolation is not forced to go through the points - this is a quality cut

        for iCoord in self.coordList:
            npoints = self.pointsArray[iCoord][:,0].shape[0]
            try:
                if npoints>0 and self.interpPresent:
                    distance = numpy.abs( self.pointsArray[iCoord][:,2] - self.doInterp(iCoord,self.pointsArray[iCoord][:,0],self.pointsArray[iCoord][:,1]) )
                    self.pointsArray[iCoord][:,3] = numpy.where(distance<maxdist,1.0,0.0)
            except:
                print("PointMesh: error in calcWgtDtoInterp")



    def calcWgtNSig(self,nsigma=3.0):
        """ set the weight for the mesh by number of sigma from TMean and TStd
        only works for interpolation type "tmean"
        wgt = 1.0 or less than nsigma, =0.0 for more...
        """

        for iCoord in self.coordList:
            npoints = self.pointsArray[iCoord][:,0].shape[0]
            try:
                if npoints>0 and self.interpTStd[iCoord]>0.0 :
                    nsigmaArray = numpy.abs(self.pointsArray[iCoord][:,2] - self.interpTMean[iCoord])/self.interpTStd[iCoord]
                    self.pointsArray[iCoord][:,3] = numpy.where(nsigmaArray<nsigma,1.0,0.0)
            except:
                print("PointMesh: error in calcWgtNSig")


    def calcWgtnMAD(self,kNN=100,nMADCut=4.0):
        """ set the weight for the mesh according to the number of MAD from the mean
        for each points.  Set the wgt = 1.0 for nMAD<cut =0.0 otherwise
        """

        listOfnMAD = []
        for iCoord in self.coordList:
            thePointsArray = self.pointsArray[iCoord]
            npoints = thePointsArray[:,0].shape[0]

            if npoints>0:

                # build the KDTree for this iCoord
                xyArr = thePointsArray[:,(0,1)]
                kdtree = cKDTree(xyArr)

                # loop over points, find kNearestNeighbors
                for ipoint in range(npoints):
                    xyPnt = thePointsArray[ipoint,(0,1)]
                    distance,neighborIndex = kdtree.query(xyPnt,kNN)

                    # passing the list of neighbors just selects those elements
                    neighborArr = thePointsArray[neighborIndex,2]

                    medianVal = numpy.median(neighborArr)
                    MAD = numpy.median(numpy.abs(neighborArr-medianVal))
                    nMAD = (thePointsArray[ipoint,2] - medianVal)/MAD
                    listOfnMAD.append(nMAD)

                    # compact pythonism to replace regular if else block
                    thePointsArray[:,3] = (1. if numpy.abs(nMAD)<nMADCut else 0.0) 
                    

        # return an array of nMAD values
        nMADArr = numpy.array(listOfnMAD)
        return nMADArr

            
    def cullMesh(self,minweight):
        """ remove points from pointsArray that have low weights
        and remake the interpolation
        """
        for iCoord in self.coordList:
            if iCoord in self.pointsArray:
                origArray = self.pointsArray[iCoord]
                npoints = origArray.shape[0]
                newList = []
                if npoints>0:
                    wgtArray = self.pointsArray[iCoord][:,3]
                    for i in range(npoints):
                        if wgtArray[i]>minweight:
                            aPoint = [origArray[i,0],origArray[i,1],origArray[i,2],origArray[i,3]]
                            newList.append(aPoint)

                newArr = numpy.array(newList)
                self.pointsArray[iCoord] = newArr.copy()

        # remake the interpolation mesh too
        if self.interpPresent:
            self.makeInterpolation(self.myMethod,self.methodVal)       


#    def cleanMesh(self,weightArray,minweight):
#        
#        # remove points from pointsArray that have low weights
#        for iCoord in self.coordList:
#            origArray = self.pointsArray[iCoord]
#            npoints = origArray.shape[0]
#            if npoints>0:
#                wgtArray = weightArray[iCoord]
#                newList = []
#                for i in range(npoints):
#                    if wgtArray[i]>minweight:
#                        aPoint = [origArray[i,0],origArray[i,1],origArray[i,2]]
#                       newList.append(aPoint)
#
#            newArr = numpy.array(newList)
#            self.pointsArray[iCoord] = newArr
#
#        # remake the interpolation mesh too
#        self.makeInterpolation(self.myMethod,self.methodVal)       


    def writePointsToFile(self,fileName):

        # convert self.pointsArray to individual numpy arrays
        sensorList = []
        xList = []
        yList = []
        zList = []
        wList = []
        for key in self.pointsArray:
            anArr = self.pointsArray[key]
            for point in anArr:
                sensorList.append(key)
                xList.append(point[0])
                yList.append(point[1])
                zList.append(point[2])
                wList.append(point[3])
        
        # output python dictionary of points to text file
        dataArray = numpy.array(list(zip(sensorList,xList,yList,zList,wList)),dtype=[('Sensor','|S3'),('x','float'),('y','float'),('z','float'),('w','float')])
        numpy.savetxt(fileName,dataArray,fmt=("%3.3s","%12.5f","%12.5f","%12.5f","%12.5f"))


    def adjustMesh(self,thetax,thetay,delta,angleconversion=1.0):
        """ adjust the values in this mesh by a rotation and offset
        """

        # loop over Coord's
        for iCoord in self.coordList:            
            # get the data
            data = self.pointsArray[iCoord]

            # make sure there is some data!
            if data.shape[0]>0 :

                # get X,Y,Z
                X = data[:,0]
                Y = data[:,1]
                Z = data[:,2]

                # remove the rotation and offset
                adjustedZ = Z - ( (X * thetay / angleconversion) + (Y * thetax / angleconversion) + delta )

                # insert this into pointsArray
                self.pointsArray[iCoord][:,2] = adjustedZ

        # remake the interpolation mesh too
        if self.interpPresent:
            self.makeInterpolation(self.myMethod,self.methodVal)       


    def mergeMesh(self,otherMesh):
        """ add another Meshes points to this one
        """
        
        # loop over Coord's
        for iCoord in self.coordList:            
            # get the data
            data = self.pointsArray[iCoord]

            # check that there is data
            if data.shape[0]>0:
                X = data[:,0]
                Y = data[:,1]
                Z = data[:,2]
                W = data[:,3]

                # get the other points
                otherData = otherMesh.pointsArray[iCoord]

                if otherData.shape[0]>0:

                    otherX = otherData[:,0]
                    otherY = otherData[:,1]
                    otherZ = otherData[:,2]
                    otherW = otherData[:,3]

                    # load new points
                    npoints = otherX.shape[0] + X.shape[0]
                    newData = numpy.zeros((npoints,4))
                    newData[:,0] = numpy.append(X,otherX)
                    newData[:,1] = numpy.append(Y,otherY)
                    newData[:,2] = numpy.append(Z,otherZ)
                    newData[:,3] = numpy.append(W,otherW)
                    self.pointsArray[iCoord] = newData.copy()

        # remake the interpolation mesh too
        if self.interpPresent:
            self.makeInterpolation(self.myMethod,self.methodVal)       

            

    def redoInterp(self,myMethod,methodVal=None):

        # reconstruct the interpolation grid for each coordinate system
        self.myMethod = myMethod
        self.methodVal = methodVal
        self.interpPresent = False
        if self.checkMethod(myMethod,methodVal):
            self.interpPresent = True
            self.makeInterpolation(self.myMethod,self.methodVal)



#    def calcDixonQ(self,inArray,outWgt):
#        """ calculate the Dixon Q ratio test for a sample of points
#        """
        
            
    def subtractMesh(self, other):
        """ subtract one mesh from another, creating a new mesh
        the new mesh uses the coordList and gridArray from
        self, but will use method='grid'
        """

        newMyMethod = 'grid'
        newCoordList = self.coordList
        newGridArray = self.gridArray
        
        # loop over gridArray points and evaluate self and other there
        newPointsDict = {}
        for iCoord in newCoordList:
            newPointsList = []
            xGrid,yGrid = self.interpGrids[iCoord]
            xxGrid = xGrid.flatten()
            yyGrid = yGrid.flatten()
            zvals = self.doInterp(iCoord,xxGrid,yyGrid) - other.doInterp(iCoord,xxGrid,yyGrid)
            for i in range(len(zvals)):
                point = [xxGrid[i],yyGrid[i],zvals[i],1.]
                newPointsList.append(point)
            newPointsArray = numpy.array(newPointsList)
            newPointsDict[iCoord] = newPointsArray    

        # build the mesh

        newMesh = PointMesh(newCoordList,newGridArray,pointsArray=newPointsDict,myMethod=newMyMethod)

        return newMesh

    def writeMesh(self,fname):
        """ write out the values at grid centers 
        """
        x,y,z,c = getGridpts()
        w = numpy.ones((z.size))
        f = open(fname)
        f.write('Sensor , x , y , z , w')
        for i in range(z.size):
            f.write("%s,%f,%f,%f,%f" % (c[i],x[i],y[i],z[i],w[i]))
        f.close()
        
    def getGridpts(self):
        """ collect all Grid points - smoothed representation of the data 
        """
        xVal = []
        yVal = []
        zVal = []
        coordVal = []
        for iCoord in self.coordList:
            xGrid,yGrid = self.interpGrids[iCoords]
            pts = self.interpValues[iCoords]
            xVal.extend(xGrid.flatten())
            yVal.extend(yGrid.flatten())
            zVal.extend(pts.flatten())
            coordVal.append(iCoord)

        xArr = numpy.array(xVal)
        yArr = numpy.array(yVal)
        zArr = numpy.array(zVal)
        coordArr = numpy.array(coordVal)
        return xArr,yArr,zArr,coordArr
        
    def calc2pt(self):
        """ calculate the two-point correlation function for a Mesh, using the values stored in interpValues
        brute force!
        """

        binSize = 10.0 # mm
        nbin = 50
        minr = 0.0
        maxr = binSize*nbin
        twoptSum = numpy.array((nbin))
        twoptNentries = numpy.array((nbin))

        x,y,z,c = self.getGridpts()
        n = x.size

        for i in range(n):
            for j in range(i,n):
                x1 = x[i]
                y1 = y[i]
                x2 = x[j]
                y2 = y[j]
                r = numpy.sqrt(numpy.power(x1-x2,2)+numpy.power(y1-y2,2))
                prod = z[i]*z[j]
                rbin = numpy.int(r/binSize)
                if rbin>=nbin:
                    twoptSum[rbin] = twoptSum[rbin] + prod
                    twoptNentries[rbin] = twoptNentires[rbin] + 1.0

        return twoptSum/twoptNentries
                    
                
        
            


        


#    def __add__(self, other):
#        """ add one mesh to another 
#        """
