#
# $Rev:: 103                                                        $:  
# $Author:: roodman                                                 $:  
# $LastChangedDate:: 2012-09-15 09:37:35 -0700 (Sat, 15 Sep 2012)   $:  
#
# ROOT Utility Library for focus and alignment code
#
import numpy
import scipy
import numpy.lib.index_tricks as itricks
import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True
import pdb
#
# declare the functions in this file
#

__all__ = ["bookandfill1d","bookfillprof","bookandfill2dfrom1dxy","bookfill2d","hshow","hfillhist","hhist","hgraph","hgraph2d","hwhisker","hcircle","hrout"]


#
# histogram book and fill routines
#
def bookandfill1d(name,title,xax,anarr):
    # book and fill a root histo
    # use 1-d Arrays xaxis,anarr
    # note that xax are bin centers

    # nb. only assumption is that xax is regularly spaced
    nx = xax.size
    xbinsiz = numpy.abs(xax[1] - xax[0])
    xmin = xax.min()-xbinsiz/2.
    xmax = xax.max()+xbinsiz/2.    
    anhist = ROOT.TH1D(name,title,nx,xmin,xmax)
    # build 1-d arrays suitable for root filling code
    anhist.FillN(nx,xax,anarr)
    return anhist

def bookfillprof(name,title,xax,anarr,n=100,xlo=0.,xhi=0.,vallo=0.,valhi=0.,opt=""):
    # book and fill a root profile histogram
    # use 1-d Arrays, or 2-d arrays, xax, anarr
    # note that xax are bin centers

    # auto find xlo,xhi 
    if xlo==0.0 and xhi==0.0 :
        xlo = xax.min()
        xhi = xax.max() + (xax.max()-xax.min())/100.

    # make the profile histogram
    if vallo==0.0 and valhi==0.0 :
        aprof = ROOT.TProfile(name,title,n,xlo,xhi,opt)
    else:
        aprof = ROOT.TProfile(name,title,n,xlo,xhi,vallo,valhi,opt)

    # if input is 2-d, convert to 1-d
    theshape = anarr.shape
    if len(theshape)==2 :
        ny,nx = theshape
        ntot = ny*nx
    elif len(theshape)==1:
        ntot = theshape[0]

    xval = xax.copy()
    yval = anarr.copy()
    xval.resize(ntot)
    yval.resize(ntot)

    # fill it
    weight = numpy.ones(ntot)
    aprof.FillN(ntot,xval,yval,weight)

    return aprof

def bookandfill2dfrom1dxy(name,title,xaxis,yaxis,anarr):
    # book and fill a root histo
    # use 1-d Arrays xaxis,yaxis,anarr
    # note that xax,yax are bin centers

    if xaxis.shape.__len__()==2 :
        xax = xaxis[0,:]
        yax = yaxis[:,0]
    else:
        xax = xaxis
        yax = yaxis
    
    nx = xax.size
    xbinsiz = numpy.abs(xax[1] - xax[0])
    xmin = xax.min()-xbinsiz/2.
    xmax = xax.max()+xbinsiz/2.    

    ny = yax.size
    ybinsiz = numpy.abs(yax[1] - yax[0])
    ymin = yax.min()-ybinsiz/2.
    ymax = yax.max()+ybinsiz/2.    

    # build 1-d arrays suitable for root filling code
    xval = xax.repeat(yax.size)
    yval = numpy.tile(yax,xax.size)

    anhist = ROOT.TH2D(name,title,nx,xmin,xmax,ny,ymin,ymax)
    nxy = nx*ny
    anarroned = anarr.copy()
    anarroned.resize(nxy)        
    #print "Input Array's shape = ", xaxis.shape, " ", yaxis.shape, " " ,anarr.shape    
    anhist.FillN(nxy,xval,yval,anarroned)
    return anhist


def bookfill2d(name,title,xaxis,yaxis,anarr):
    # book and fill a root histo
    # use 2-d Arrays xaxis,yaxis,anarr only

    xax = xaxis[0,:]
    yax = yaxis[:,0]
    
    nx = xax.size
    xbinsiz = numpy.abs(xax[1] - xax[0])
    xmin = xax.min()-xbinsiz/2.
    xmax = xax.max()+xbinsiz/2.    

    ny = yax.size
    ybinsiz = numpy.abs(yax[1] - yax[0])
    ymin = yax.min()-ybinsiz/2.
    ymax = yax.max()+ybinsiz/2.    

    anhist = ROOT.TH2D(name,title,nx,xmin,xmax,ny,ymin,ymax)
    nxy = nx*ny
    anarroned = anarr.copy()
    xaxised = xaxis.copy()
    yaxised = yaxis.copy()
    anarroned.resize(nxy)        
    xaxised.resize(nxy)
    yaxised.resize(nxy)
    #print "Input Array's shape = ", xaxis.shape, " ", yaxis.shape, " " ,anarr.shape    
    anhist.FillN(nxy,xaxised,yaxised,anarroned)
    return anhist


def hshow(name,title,anarr,opt="zcol"):
    # book and fill a root histo
    # use 1-d and 2-d Arrays anarr

    # fancy trick to keep histos and canvas in a permanent scope!
    global histlist
    global canvaslist
    try:
        if histlist is None:
            histlist = []
    except NameError:
        histlist = []
    try:
        if canvaslist is None:
            canvaslist = []
    except NameError:
        canvaslist = []

    # decide on 1-d or 2-d here
    shape = anarr.shape

    if len(shape)==2 :
        ny,nx = anarr.shape
        yaxis,xaxis = itricks.mgrid[0.0:float(ny):1.,0.0:float(nx):1.]

        xax = xaxis[0,:]
        yax = yaxis[:,0]

        nx = xax.size
        xbinsiz = numpy.abs(xax[1] - xax[0])
        xmin = xax.min()-xbinsiz/2.
        xmax = xax.max()+xbinsiz/2.    

        ny = yax.size
        ybinsiz = numpy.abs(yax[1] - yax[0])
        ymin = yax.min()-ybinsiz/2.
        ymax = yax.max()+ybinsiz/2.    

        anhist = ROOT.TH2D(name,title,nx,xmin,xmax,ny,ymin,ymax)
        nxy = nx*ny

        anarroned = numpy.float64(anarr.copy())
        xaxised = xaxis.copy()
        yaxised = yaxis.copy()
        anarroned.resize(nxy)        
        xaxised.resize(nxy)
        yaxised.resize(nxy)
        #print "Input Array's shape = ", xaxis.shape, " ", yaxis.shape, " " ,anarr.shape
        anhist.FillN(nxy,xaxised,yaxised,anarroned)
        cName = "c" + name
        c = ROOT.TCanvas(cName,cName)
        anhist.Draw(opt)
        histlist.append(anhist)
        canvaslist.append(c)

    elif len(shape)==1 :

        nx = anarr.shape[0]
        xax = itricks.mgrid[0.0:float(nx):1.]
        
        nx = xax.size
        xbinsiz = numpy.abs(xax[1] - xax[0])
        xmin = xax.min()-xbinsiz/2.
        xmax = xax.max()+xbinsiz/2.    

        anhist = ROOT.TH1D(name,title,nx,xmin,xmax)
        anarred = anarr.copy()   #FillN doesn't always work, unless we have local copy
        anhist.FillN(nx,xax,anarred,1)
        cName = "c" + name
        c = ROOT.TCanvas(cName,cName)
        anhist.Draw()
        histlist.append(anhist)
        canvaslist.append(c)
        
    
    return anhist,c


def hhist(name,title,anarr,nbins=100,lo=-999,hi=-999):
    # book and fill a root histo
    # use 1-d and 2-d Arrays anarr as input
    # histogram the vlaues

    # fancy trick to keep histos and canvas in a permanent scope!
    global histlist
    global canvaslist
    try:
        if histlist is None:
            histlist = []
    except NameError:
        histlist = []
    try:
        if canvaslist is None:
            canvaslist = []
    except NameError:
        canvaslist = []

    # put in 1d array
    shape = anarr.shape
    if len(shape)==2 :
        ny,nx = anarr.shape
        nxy = nx*ny
        anarroned = anarr.copy()
        anarroned.resize(nxy)
    else :
        anarroned = anarr.copy()
        nxy = anarr.shape[0]

    # book it
    if lo==hi :
        lo = numpy.min(anarroned)*0.95
        hi = numpy.max(anarroned)*1.05
    anhist = ROOT.TH1D(name,title,nbins,lo,hi)

    # fill it
    wgt = numpy.ones(nxy)
#    anhist.FillN(nxy,anarroned,wgt)
    for i in range(nxy):
        anhist.Fill(anarroned[i])
        
    cName = "c" + name
    c = ROOT.TCanvas(cName,cName)
    anhist.Draw()
    histlist.append(anhist)
    canvaslist.append(c)
        
    return anhist,c


def hfillhist(name,title,anarr,nbins=100,lo=-999,hi=-999):
    # book and fill a root histo
    # use 1-d and 2-d Arrays anarr as input
    # histogram the vlaues

    # put in 1d array
    shape = anarr.shape
    if len(shape)==2 :
        ny,nx = anarr.shape
        nxy = nx*ny
        anarroned = anarr.copy()
        anarroned.resize(nxy)
    else :
        anarroned = anarr.copy()
        nxy = anarr.shape[0]

    # book it
    if lo==hi :
        lo = numpy.min(anarroned)*0.95
        hi = numpy.max(anarroned)*1.05
    anhist = ROOT.TH1D(name,title,nbins,lo,hi)

    # fill it
    wgt = numpy.ones(nxy)
    anhist.FillN(nxy,anarroned,wgt)
    return anhist


def hgraph(name,title,xarr,yarr,eyarr=None,opt="AP"):
    # load a TGraph

    # fancy trick to keep histos and canvas in a permanent scope!
    global histlist
    global canvaslist
    try:
        if histlist is None:
            histlist = []
    except NameError:
        histlist = []
    try:
        if canvaslist is None:
            canvaslist = []
    except NameError:
        canvaslist = []

    if eyarr!=None :
        exarr = numpy.zeros((xarr.size))
        anhist = ROOT.TGraphErrors(xarr.size,xarr,yarr,exarr,eyarr)
    else:
        anhist = ROOT.TGraph(xarr.size,xarr,yarr)

    anhist.SetTitle(title)
    anhist.SetName(name)

    cName = "c" + name
    c = ROOT.TCanvas(cName,cName)
    anhist.Draw(opt)

    histlist.append(anhist)
    canvaslist.append(c)
        
    return anhist,c


def hgraph2d(name,title,xarr,yarr,zarr,opt="zcolpcol"):
    # load a TGraph2D

    # fancy trick to keep histos and canvas in a permanent scope!
    global histlist
    global canvaslist
    try:
        if histlist is None:
            histlist = []
    except NameError:
        histlist = []
    try:
        if canvaslist is None:
            canvaslist = []
    except NameError:
        canvaslist = []

    anhist = ROOT.TGraph2D(xarr.size,xarr,yarr,zarr)

    anhist.SetTitle(title)
    anhist.SetName(name)

    cName = "c" + name
    c = ROOT.TCanvas(cName,cName)
    anhist.Draw(opt)

    histlist.append(anhist)
    canvaslist.append(c)
        
    return anhist,c

def hwhisker(name,title,xmin,xmax,ymin,ymax,whiskermax,scale,npt,xval,yval,magval,angleval):
    # make a whisker plot, given arrays with x,y,magnitude,angle for discrete points
    
    # fancy trick to keep histos and canvas in a permanent scope!
    global histlist
    global canvaslist
    global linelist
    try:
        if histlist is None:
            histlist = []
    except NameError:
        histlist = []
    try:
        if canvaslist is None:
            canvaslist = []
    except NameError:
        canvaslist = []
    try:
        if linelist is None:
            linelist = []
    except NameError:
        linelist = []

    cName = "c" + name
    c = ROOT.TCanvas(cName,cName)
    canvaslist.append(c)

    hframe = ROOT.TH2D("hframe","",100,xmin,xmax,100,ymin,ymax);
    histlist.append(hframe)    
    hframe.Draw();
    hframe.SetStats(0)
        
    # calculate whisker x,y lengths
    whiskerYsiz = scale* magval * numpy.sin(angleval) / whiskermax
    whiskerXsiz = scale* magval * numpy.cos(angleval) / whiskermax

    whiskerYlo = yval - 0.5*whiskerYsiz
    whiskerYhi = yval + 0.5*whiskerYsiz
    whiskerXlo = xval - 0.5*whiskerXsiz
    whiskerXhi = xval + 0.5*whiskerXsiz

    
    #draw whiskers
    linelist = []
    for ipt in range(npt):
        
        linelist.append(ROOT.TLine(whiskerXlo[ipt],whiskerYlo[ipt],whiskerXhi[ipt],whiskerYhi[ipt]))
        linelist[ipt].Draw("same")

    c.Update()

    return c

    

def hcircle(x,y,r):
    # load a TEllipse

    # fancy trick to keep objects in a permanent scope!
    global histlist
    try:
        if histlist is None:
            histlist = []
    except NameError:
        histlist = []

    tcc = ROOT.TEllipse(x,y,r)
    tcc.SetFillStyle(1)
    tcc.SetLineColor(2)
    tcc.SetLineWidth(3)
    tcc.Draw()

    histlist.append(tcc)        
    return tcc

def hrout(fileName):

    currentDir = ROOT.gDirectory
    myList = ROOT.gDirectory.GetList()
    output = ROOT.TFile(fileName,"RECREATE")
    output.cd()
    for obj in myList:
        if not obj.InheritsFrom("TTree")  :
            print(obj.GetName(), obj.GetTitle())
            obj.Write()
                
    output.Close()
    currentDir.cd()




