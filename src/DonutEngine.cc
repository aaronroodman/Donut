//
// DonutEngine.cc:  Engine to calculate focal plane image from 
//                  pupil plane Zernike expansion
// Copyright (C) 2011 Aaron J. Roodman, SLAC National Accelerator Laboratory, Stanford University
//

// C and C++ headers
#include <iostream>
#include <cmath>
#include <string>
#include <time.h>

// Class header files
#include "fitsio.h"
#include "DonutEngine.h"
#include "FFTWClass.h"

// Constructors

DonutEngine::DonutEngine(MapStoS inputMapS, MapStoI inputMapI, MapStoD inputMapD){

  // save these input options
  _inputMapS = inputMapS;
  _inputMapI = inputMapI;
  _inputMapD = inputMapD;
  
  fillOptions(inputMapS,inputMapI,inputMapD);

  // print if desired
  printOptions();

  // all initialization
  init();

}

//DonutEngine::DonutEngine(int iTelescope, 
//			 double waveLength, 
//			 int nZernikeTerms, 
//			 int nbin, 
//			 int nPixels, 
//			 const char* outputPrefix,
//			 bool debugFlag,
//			 int printLevel,
//			 bool gridCalcMode,
//			 int pixelOverSample,
//			 double scaleFactor,
//			 const char* inputPupilMask,
//			 int zemaxToDECamSignFlip,
//			 bool calcRzeroDerivative){

//   _iTelescope = iTelescope;
//   _waveLength = waveLength;
//   _nZernikeTerms = nZernikeTerms;  //default was 37
//   _nbin = nbin; // scale if desired
//   _nPixels = nPixels;
//   _outputPrefix = outputPrefix;
//   _debugFlag = debugFlag;
//   _printLevel = printLevel;
//   _gridCalcMode = gridCalcMode;
//   _pixelOverSample = pixelOverSample;   // scale if desired
//   _scaleFactor = scaleFactor;  // scale if desired
//   _inputPupilMask = inputPupilMask;
//   _zemaxToDECamSignFlip = zemaxToDECamSignFlip;  // use =-1 for x(Zemax) => x(DECam FITs)  AND y(Zemax)=> y(DECam FITs) 
//   _calcRzeroDerivative = calcRzeroDerivative;
//   //   added y flip in code below on 10/3/2012, AJR

//   // always initialize xDECam,yDECam to zero, change with setXYDECam
//   _xDECam = 0.0;
//   _yDECam = 0.0;

//   init();

//   printOptions();
// }

DonutEngine::~DonutEngine(){

  if (_debugFlag){
    int status = 0;         /* initialize status before calling fitsio routines */
    if (fits_close_file(_fptr, &status)){            /* close the file */
      fits_report_error(stderr, status);
    }
  }
  delete _fft2PlanC;
  delete _ifft2PlanC;
}

void DonutEngine::closeFits(){

  if (_debugFlag){
    int status = 0;         /* initialize status before calling fitsio routines */
    if (fits_close_file(_fptr, &status)){            /* close the file */
      fits_report_error(stderr, status);
    }
  }

}

void DonutEngine::fillOptions(MapStoS inputMapS, MapStoI inputMapI, MapStoD inputMapD){

  // default values for input maps
  MapStoS defaultMapS;
  defaultMapS["outputPrefix"] = "test";
  defaultMapS["inputPupilMask"] = "";

  MapStoI defaultMapI;
  defaultMapI["iTelescope"] = 0;
  defaultMapI["nbin"] = 256;
  defaultMapI["nPixels"] = 64;
  defaultMapI["pixelOverSample"] = 4;
  defaultMapI["nZernikeTerms"] = 11;
  defaultMapI["printLevel"] = 0;
  defaultMapI["debugFlag"] = 0;
  defaultMapI["gridCalcMode"] = 1;
  defaultMapI["zemaxToDECamSignFlip"] = 1;   //CHANGED DEFAULT to positive 1 on 10/4/2012 AJR
  defaultMapI["calcRzeroDerivative"] =0;

  MapStoD defaultMapD;
  defaultMapD["waveLength"] = 700.0e-9;
  defaultMapD["scaleFactor"] = 2.0;

  // loop over maps and insert input values
  MapStoS optionMapS;
  for (MapStoS::iterator i = defaultMapS.begin(); i != defaultMapS.end(); ++i){
    std::string key = i->first;
    // is this key in the inputMap?
    if (inputMapS.count(key)==0){
      optionMapS[key] = defaultMapS[key];
    } else {
      optionMapS[key] = inputMapS[key];
    }
  }

  MapStoI optionMapI;
  for (MapStoI::iterator i = defaultMapI.begin(); i != defaultMapI.end(); ++i){
    std::string key = i->first;
    // is this key in the inputMap? check both inputMapI and inputMapD
    if (inputMapI.count(key)!=0){
      optionMapI[key] = inputMapI[key];
    } else if (inputMapD.count(key)!=0) {
      optionMapI[key] = int(inputMapD[key]+0.5);  //in case we specify an integer using a float
    } else {
      optionMapI[key] = defaultMapI[key];
    }
  }

  MapStoD optionMapD;
  for (MapStoD::iterator i = defaultMapD.begin(); i != defaultMapD.end(); ++i){
    std::string key = i->first;
    // is this key in the inputMap? check both inputMapI and inputMapD
    if (inputMapD.count(key)!=0){
      optionMapD[key] = inputMapD[key];
    } else if (inputMapI.count(key)!=0) {
      optionMapD[key] = inputMapI[key];  //in case we specify a float using an integer
    } else {
      optionMapD[key] = defaultMapD[key];
    }
  }


  // fill variables from combined maps
  _iTelescope = optionMapI["iTelescope"];
  _nZernikeTerms = optionMapI["nZernikeTerms"];
  _nbin = optionMapI["nbin"];
  _nPixels = optionMapI["nPixels"];
  _debugFlag = bool(optionMapI["debugFlag"]);
  _printLevel = optionMapI["printLevel"];
  _pixelOverSample = optionMapI["pixelOverSample"];
  _gridCalcMode = bool(optionMapI["gridCalcMode"]);
  _zemaxToDECamSignFlip = optionMapI["zemaxToDECamSignFlip"];
  _calcRzeroDerivative = optionMapI["calcRzeroDerivative"];

  _waveLength = optionMapD["waveLength"];
  _scaleFactor = optionMapD["scaleFactor"];

  _outputPrefix = optionMapS["outputPrefix"];
  _inputPupilMask = optionMapS["inputPupilMask"];

  // always initialize xDECam,yDECam to zero, change with setXYDECam
  _xDECam = 0.0;
  _yDECam = 0.0;
  //_xDESI  = 0.0;
  //_yDESI  = 0.0;

}

void DonutEngine::printOptions(){

  if (_printLevel>=2){
    std::cout << "iTelescope = " << _iTelescope << std::endl;
    std::cout << "waveLength = " << _waveLength << std::endl;
    std::cout << "nZernikeTerms = " << _nZernikeTerms << std::endl;
    std::cout << "nbin = " << _nbin << std::endl;
    std::cout << "nPixels = " << _nPixels << std::endl;
    std::cout << "outputPrefix = " << _outputPrefix << std::endl;
    std::cout << "debugFlag = " << _debugFlag << std::endl;
    std::cout << "printLevel = " << _printLevel << std::endl;
    std::cout << "gridCalcMode = " << _gridCalcMode << std::endl;
    std::cout << "pixelOverSample = " << _pixelOverSample << std::endl;
    std::cout << "scaleFactor = " << _scaleFactor << std::endl;
    std::cout << "inputPupilMask = " << _inputPupilMask << std::endl; 
    std::cout << "zemaxToDECamSignFlip = " << _zemaxToDECamSignFlip << std::endl; 
    std::cout << "calcRzeroDerivative = " << _calcRzeroDerivative << std::endl; 
  }


}

void DonutEngine::init(){
  
  _M_PI = 3.1415926535897931;  //replace with math.h ??
  _M_PI_4 = _M_PI/4.;

  // used for alignment of arrays
  alignR=sizeof(Real);
  alignC=sizeof(Complex);

  // if debug is on, open a fits file
  if (_debugFlag){
    //    const char* filename = {"testcpluplus.fits"};
    std::string filename = std::string(_outputPrefix) + "-debug.fits";
    remove(filename.c_str());
    int status = 0;         /* initialize status before calling fitsio routines */
    if (fits_create_file(&_fptr,filename.c_str(), &status)){   /* create new file */
      fits_report_error(stderr, status);
    }
  }
    
  // setup arrays, parameters for DonutEngine  
  calcParameters(_iTelescope);
  setupArrays();
  setupStuff();
  initStateMachine();
  defineParams();

  // internal diagnostics
  resetTimers();
  nCallsCalcAll = 0;
  nCallsCalcDerivative = 0;

  // we may want to calculate the Rzero Derivative, so make another DonutEngine to aid this calculation
  if (_calcRzeroDerivative) {
    MapStoS inS = _inputMapS;
    MapStoI inI = _inputMapI;
    MapStoD inD = _inputMapD;

    inS["outputPrefix"] = "calcDerivative";
    inI["calcRzeroDerivative"] = 0;  //Turn this off, avoid an infinite loop!
     
    _anotherDonutEngine  = new DonutEngine(inS,inI,inD);
  }

}

void DonutEngine::calcParameters(int iT){

  // physical parameters
  // iTelescope = 0 for Blanco + DECam
  //            = 1 for Blanco + MosaicII
  //            = 2 for LSST
  //            = 3 for Magellan MegaCam
  //            = 4 for Magellan IMACS F/2
  //            = 5 for DESI


  if (iT==0) {
    _outerRadius = 0.7174;  // see AJR logbook, DES Vol3, pg 25
    _innerRadius = 0.301;
    _zLength = 4.274419;
    _lambdaz = _waveLength * _zLength;  //see DES notebook Vol1 pg 78-80
    _fLength = 11.719;
    _pixelSize = 15.0e-6;
  } else if (iT==1){   
    _outerRadius = 0.29158;
    _innerRadius = 0.122;
    _zLength = 1.696231;
    _lambdaz = _waveLength * _zLength;  
    _fLength = 11.45784;
    _pixelSize = 15.0e-6; 
  } else if (iT==2){
    _outerRadius = 1.10874;   // this the radius of the pupil as found in Zemax!
    _innerRadius = 0.61 * _outerRadius;   // 4/4/2016 look at Zemax gives 0.612 - use 0.61 for LSST donuts
    _zLength = 2.734706;
    _lambdaz = _waveLength * _zLength;  
    _fLength = 10.31007;
    _pixelSize = 10.0e-6;
  } else if (iT==3){
    _outerRadius = 6.5/2.0;
    _innerRadius = 0.352 * _outerRadius; // from Povilas
    _zLength = 34.97;  // F#5.38 from Povilas
    _lambdaz = _waveLength * _zLength;  
    _fLength = 34.97;
    _pixelSize = 2.0 * 13.5e-6;  // operated in 2x2 mode
  } else if (iT==4){
    _outerRadius = 6.5/2.0;
    _innerRadius = 0.352 * _outerRadius; // from Povilas
    _zLength = 15.47;  // F#2.38 see http://www.lco.cl/telescopes-information/magellan/instruments/imacs/imacs-specs
    _lambdaz = _waveLength * _zLength;  
    _fLength = 15.47;
    _pixelSize = 1.0 * 15.0e-6;  // operated in 1x1 mode, 15 micron pixels
  } else if (iT==5){
    _outerRadius = 0.892; // from Zemax
    _innerRadius = 0.439; // from outerRadius*<innermajoraxis, innerminoraxis>/<outermajoraxis, outerminoraxis>
    _zLength = 6.523858;  // F#3.66
    _lambdaz = _waveLength * _zLength;  
    _fLength = 13.92452; //from Zemax
    _pixelSize = 1.0 * 15.0e-6;  // operated in 1x1 mode, 15 micron pixels

  } 

  //
  //  here we fix the size of the grid bins on the focal plane to be equal to the 
  //  pixelsize/pixelOverSample. 
  //   (see pg 131 DES vol1 and FFTnotes.pdf)
  //
  // NOTE added 10/3/2012:  Code will break if _scaleFactor is < 1.0
  // need to manually adjust scaleFactor to ensure that _scaleFactor > 1.0

  Real F =  _zLength/(2.* _outerRadius);
  Real lambdaF =  F * _waveLength;
  
  if (_gridCalcMode){    
    _ngridperPixel =  _pixelOverSample;
    _deltaX =  _pixelSize / _ngridperPixel;
    _pupilscale = ( lambdaF / _deltaX) *  _scaleFactor;
    _Lu =  2.0 *  _outerRadius * _pupilscale;
  } else {
    // when we start from WFM->Focal plane, we need to specify Lu first, 
    // then find deltaX from it
    // in this case we do interpolation instead of just using every nth grid point
    _Lu = 2.0 * _outerRadius;
    _deltaX =  lambdaF;
    _ngridperPixel =  _pixelSize / _deltaX;
  }

  //
  // other needed things
  //
  _Lf = 1./ _deltaX;
  _nhalfPixels =  _nPixels/2;

  // print output
  if (_printLevel>=2){
    std::cout << "ngridperPixel = " << _ngridperPixel << std::endl;
    std::cout << "deltaX = " << _deltaX << std::endl;
    std::cout << "pupilscale = " << _pupilscale << std::endl;
    std::cout << "Lu = " << _Lu << std::endl;
    std::cout << "Lf = " << _Lf << std::endl;
    std::cout << "nhalfPixels = " << _nhalfPixels << std::endl;
  }

}

void DonutEngine::itricksMGrid(int n,Real lo, Real hi, Matrix& x, Matrix& y){
  // form an x,y grid to mimic python itricks.mgrid
  // assume that it is a regular square grid

  x.Dimension(n,n);
  y.Dimension(n,n);
  x.Activate(alignR);
  y.Activate(alignR);

  Real delta = (hi-lo)/((Real)n - 1. );
  for (int j=0;j<n;j++){
    for (int i=0;i<n;i++){

      //row major order offset = row*NUMCOLS + column
      // but arrays are Arr[iy][ix]

      x(j,i) = lo + (Real) i * delta;
      y(j,i) = lo + (Real) j * delta;

    }
  }

}

void DonutEngine::makePupilArrays(int nbin,Real lo,Real hi,Real radiusOuter){
  // make the rho,theta,xaxis,yaxis,pupilfunc arrays here

  // build the x,y grid, defined so that imagearray[Y,X]
  // symmetric around 0.0 now!!!
  itricksMGrid(nbin,lo,hi,_xaxis,_yaxis);

  // cartesian -> polar
  _rho.Dimension(nbin,nbin);
  _theta.Dimension(nbin,nbin);

  _rho.Activate(alignR);
  _theta.Activate(alignR);

  for (int i=0;i<nbin*nbin;i++){
    _rho(i) = sqrt(_xaxis(i)*_xaxis(i) + _yaxis(i)*_yaxis(i))/radiusOuter;  
    _theta(i) = atan2(_yaxis(i),_xaxis(i));
  }

  
}

void DonutEngine::makeXPsf(int nbin,Real lambdaz){
  // build grid for focal plane, from grid for pupil plane

  // see DES logbook Vol1, page 25
  // and page 87 and 90

  //  Real deltax = _xaxis(0,1) - _xaxis(0,0);
  Real deltay = _yaxis(1,0) - _yaxis(0,0);
  // ASSUME deltax==deltay !!!

  // now symmetric around 0.0 (check this with logbook)
  itricksMGrid(nbin,-1./(2.*deltay),1./(2.*deltay),_xpsf,_ypsf);
  _xpsf *= lambdaz;
  _ypsf *= lambdaz;
    
}

void DonutEngine::setupArrays(){

  // pupil Mask
  _pupilMask.Dimension(_nbin,_nbin);
  _pupilMask.Activate(alignR);
  _pupilSNorm = 0.0;

  // pupilWave array
  _pupilWaveZernike.Dimension(_nbin,_nbin);
  _pupilWaveZernike.Activate(alignR);

  // MUST explicitly zero _pupilWaveZernike 
  for (int i=0;i<_nbin*_nbin;i++){
    _pupilWaveZernike(i) = 0.0;
  }

  // pupilWave array
  _pupilWaveZernikePlusDelta.Dimension(_nbin,_nbin);
  _pupilWaveZernikePlusDelta.Activate(alignR);

  // MUST explicitly zero _pupilWaveZernikePlusDelta
  for (int i=0;i<_nbin*_nbin;i++){
    _pupilWaveZernikePlusDelta(i) = 0.0;
  }

  // set dimensionality of pupil arrays
  _pupilFunc.Dimension(_nbin,_nbin);
  _pupilFunc.Activate(alignC);
  _pupilFuncStar.Dimension(_nbin,_nbin);
  _pupilFuncStar.Activate(alignC);

  // Psf array
  _calcG.Dimension(_nbin,_nbin);
  _calcG.Activate(alignC);
  _calcGstar.Dimension(_nbin,_nbin);
  _calcGstar.Activate(alignC);
  _magG.Dimension(_nbin,_nbin);
  _magG.Activate(alignC);
  _phaseG.Dimension(_nbin,_nbin);
  _phaseG.Activate(alignC);
  _psfOptics.Dimension(_nbin,_nbin);
  _psfOptics.Activate(alignR);
  _ftsOptics.Dimension(_nbin,_nbin);
  _ftsOptics.Activate(alignC);

  // Atmos
  _rAtmos.Dimension(_nbin,_nbin);
  _rAtmos.Activate(alignR);
  _shftrAtmos.Dimension(_nbin,_nbin);
  _shftrAtmos.Activate(alignR);
  _psfAtmos.Dimension(_nbin,_nbin);
  _psfAtmos.Activate(alignR);
  _ftsAtmos.Dimension(_nbin,_nbin);
  _ftsAtmos.Activate(alignC);

  // arrays for Pixelization
  _ftsPixel.Dimension(_nbin,_nbin);
  _ftsPixel.Activate(alignC);

  // arrays for convolution
  _convOpticsAtmosPixel.Dimension(_nbin,_nbin);
  _convOpticsAtmosPixel.Activate(alignR);

  // arrays for Pixel datq
  _valPixelCenters.Dimension(_nPixels,_nPixels);
  _valPixelCenters.Activate(alignR);
  _calcImage.Dimension(_nPixels,_nPixels);
  _calcImage.Activate(alignR);

  // arrays for deltaWFM
  _deltaWFM.Dimension(_nbin,_nbin);
  _deltaWFM.Activate(alignR);
    
}


void DonutEngine::setupStuff(){

  // Fourier Transform arrays

  _fftInputArray.Dimension(_nbin,_nbin);
  _fftOutputArray.Dimension(_nbin,_nbin);
  _ifftInputArray.Dimension(_nbin,_nbin);
  _ifftOutputArray.Dimension(_nbin,_nbin);
  _fftrtcInputArray.Dimension(_nbin,_nbin);
  _fftrtcTempArray.Dimension(_nbin,(_nbin/2) + 1);
  _fftrtcOutputArray.Dimension(_nbin,_nbin);

  _fftInputArray.Activate(alignC);
  _fftOutputArray.Activate(alignC);
  _ifftInputArray.Activate(alignC);
  _ifftOutputArray.Activate(alignC);
  _fftrtcInputArray.Activate(alignR);
  _fftrtcTempArray.Activate(alignC);
  _fftrtcOutputArray.Activate(alignC);

  

  // need for debugging
  //_fftrtcInputArray.AtIndex(1,1);
  //_fftrtcOutputArray.AtIndex(1,1);

  _fft2PlanC =  new fftw2dctc(_fftInputArray,_fftOutputArray,-1);
  _ifft2PlanC = new fftw2dctc(_ifftInputArray,_ifftOutputArray,1);
  _fft2rtcPlanC = new fftw2drtc(_fftrtcInputArray,_fftrtcTempArray,_fftrtcOutputArray);

  // setup the Pupil function and PSF arrays here
  makePupilArrays(_nbin,-_Lu/2.0,_Lu/2.0,_outerRadius);
  makeXPsf(_nbin,_scaleFactor*_lambdaz);

  // make Zernike basis
  _zernikeObject = new Zernike(_rho,_theta,_nZernikeTerms);

  // FT of Pixel-sized box - just need this once
  calcFTPixel();

  // arrays for Atmosphere
  Matrix xAtmos,yAtmos;
  itricksMGrid(_nbin,-_Lf/2.,_Lf/2.,xAtmos,yAtmos);

  for (int i=0;i<_nbin*_nbin;i++){
    _rAtmos(i) = sqrt(xAtmos(i)*xAtmos(i) + yAtmos(i)*yAtmos(i));
  }
  _shftrAtmos = _rAtmos;    // deep copy 
  fftShift(_shftrAtmos); // shifts in place, was InvShift

  // initialize flag for deltaWFM
  _usingDeltaWFM = false;
  
}

void DonutEngine::initStateMachine(){
        
  // State machine flags
  _first = true;
  _statePupilMask = false;
  _statePupilFunc = false;
  _stateAtmos = false;
  _stateDerivatives = false;

}


void DonutEngine::defineParams(){
  // define parameters here

  // always omit the piston term - so if nZernikeTerms=37, we only use 36
  nZernikeSize = _nZernikeTerms - 1;

  // build array for Zernike parameters
  _ZernikeArr.Dimension(nZernikeSize);
  _ZernikeArr.Activate();

  // build pointers to parameters here
  int ipar = -1;
  ipar_nEle     = ipar = ipar + 1;  
  ipar_rzero	 = ipar = ipar + 1;
  ipar_bkgd     = ipar = ipar + 1;
  ipar_ZernikeFirst = ipar = ipar + 1;
  ipar_ZernikeLast  = ipar = ipar + nZernikeSize - 1;
  npar = ipar_ZernikeLast+1;

  // define parameter names here
  parNames.resize(npar);
  parTitles.resize(npar);
    
  parNames[ipar_nEle] = "nele";
  parNames[ipar_rzero] = "rzero";
  parNames[ipar_bkgd] = "bkgd";
  parTitles[ipar_nEle] = "Number of photo-electrons";
  parTitles[ipar_rzero] = "Fried Parameter [m]";
  parTitles[ipar_bkgd] = "Background [npe]";

  for (int iZ=0;iZ<nZernikeSize;iZ++){
    std::ostringstream aName;
    aName << "zern" << iZ+2;
    parNames[ipar_ZernikeFirst+iZ] = aName.str();
    if ((iZ+1)<37){
      parTitles[ipar_ZernikeFirst+iZ] = _zernikeObject->_zernikeDescription[iZ+1];
    }
  }

  // array for current parameters 
  _parCurrent.Dimension(npar);
  _parCurrent.Activate();

  // Derivatives
  _dChi2dpar.Dimension(npar);
  _dChi2dpar.Activate();

  // zero them
  for (int ipar=0;ipar<npar;ipar++){
    _parCurrent[ipar] = 0.0;
    _dChi2dpar[ipar] = 0.0;
  }

}

void DonutEngine::calcWFMtoImage(double* wfm, int nx, int ny){
  Matrix wfmM(nx,ny);
  for (int i=0;i<nx*ny;i++){
    wfmM(i) = wfm[i];
  }
  calcWFMtoImage(wfmM);
}
 
void DonutEngine::calcWFMtoImage(Matrix& wfm){
      
  // need to fill Par separately

  _statePupilFunc = false; 
      
  // State Machine: call each step of the calculation
  if (!_statePupilMask){ 
    calcPupilMask();
  } 
  if ( (!_statePupilMask) || (!_statePupilFunc)){
    calcPupilFuncFromWFM(wfm);
    calcOptics();
  } 
  if (!_stateAtmos){ 
    calcAtmos();
  } 
  if ( (!_statePupilMask) || (!_statePupilFunc) || (!_stateAtmos) ){
    calcConvolute();
    calcPixelate();
  }

  // calculate image, also in electrons
  for (int i=0;i<_nPixels*_nPixels;i++){
    _calcImage(i) = _nEle*_valPixelCenters(i) + _bkgd;
  }

}

void DonutEngine::setDeltaWFM(double* wfm, int nx, int ny){
  Matrix wfmM(nx,ny);
  for (int i=0;i<nx*ny;i++){
    wfmM(i) = wfm[i];
  }
  setDeltaWFM(wfmM);
}
  
void DonutEngine::setDeltaWFM(Matrix& wfm){
      
  // set the delta WFM contribution

  // this flag is used in fillPars and calcPupilFunction
  _usingDeltaWFM = true;

  // explictly set state flag
  _statePupilFunc = false;

  // fill
  for (int i=0;i<_nbin*_nbin;i++){
    _deltaWFM(i) = wfm(i);
  }

}

void DonutEngine::unsetDeltaWFM(){
  _usingDeltaWFM = false;

  // explictly set state flag
  _statePupilFunc = false;

    // fill
  for (int i=0;i<_nbin*_nbin;i++){
    _deltaWFM(i) = 0.0;
  }

}


void DonutEngine::calcPupilFuncFromWFM(Matrix& wfm){

  clock_t start = clock();

  if (_printLevel>=2){
    std::cout << "DonutEngine: calcPupilFuncFromWFM " << std::endl;
  }
  
  // calculate the pupilFunc(tion) from the WFM
  //  Complex I = Complex(0.0,1.0);
  Complex twopiI = Complex(0.0,2.0*_M_PI);
  for (int i=0;i<_nbin*_nbin;i++){
    if (_pupilMask(i)==0.0){
      _pupilWaveZernike(i) = 0.0;
      _pupilFunc(i) = 0.0;
      _pupilFuncStar(i) = 0.0;
    } else {
      _pupilWaveZernike(i) = wfm(i);
      _pupilFunc(i) = _pupilMask(i) * exp(twopiI  * wfm(i));   // no lambda here, so units are in waveLength
      _pupilFuncStar(i) = conj(_pupilFunc(i));
    }
  }

  if (_debugFlag && nCallsCalcAll<=1){    
    toFits(_fptr,wfm);
    toFits(_fptr,_pupilFunc);
  }

  clock_t stop = clock();
  _timePupilFunc += (stop-start)/(Real)CLOCKS_PER_SEC;

}


void DonutEngine::calcAll(double* par, int n){  
  calcAll(par);
}

void DonutEngine::calcAll(Real* par){
  
  nCallsCalcAll++;

  // fill parameters and determine how much of the calculation to repeat 
  fillPar(par);
        
  // State Machine: call each step of the calculation
  if (!_statePupilMask){ 
    calcPupilMask();
  } 
  if ( (!_statePupilMask) || (!_statePupilFunc) ){ 
    calcPupilFunc();
    calcOptics();
  } 
  if (!_stateAtmos){ 
    calcAtmos();
  } 
  if ( (!_statePupilMask) || (!_statePupilFunc) || (!_stateAtmos) ){ 
    calcConvolute();
    calcPixelate();
  }

  // calculate image
  for (int i=0;i<_nPixels*_nPixels;i++){
    _calcImage(i) = _nEle * _valPixelCenters(i) + _bkgd;
  }

  // calculate Derivatives
  //        calcDerivatives(image,weight)

  // save the parameters
  savePar();
            
  if (_printLevel>=2){
    std::cout << "DonutEngine: calcAll is done" << std::endl;
  }

}


void DonutEngine::calcAll(double* par, int n, double* dwfm, int nx, int ny){
  Matrix dwfmM(nx,ny);
  for (int i=0;i<nx*ny;i++){
    dwfmM(i) = dwfm[i];
  }

  // fill deltaWFM  - need to do this before fillPar inside calcAll
  setDeltaWFM(dwfmM);

  // now just call calcAll
  calcAll(par);
}


void DonutEngine::fillPar(double* par, int n){
  fillPar(par);
}

void DonutEngine::fillPar(Real* par){
  // save input parameters, and also reevaluate State Machine based on values of par
  // and on values stored in _xDECam,_yDECam


  // save input parameters
  for (int ipar=0;ipar<npar;ipar++){
    _parCurrent[ipar] = par[ipar];
  }

  // store fit parameters
  _nEle = par[ipar_nEle];
  //
  // scale Zernike array by scaleFactor
  //
  _rzero = par[ipar_rzero];
  _bkgd = par[ipar_bkgd];
  for (int iZ=0;iZ<nZernikeSize;iZ++){
    _ZernikeArr[iZ] = par[ipar_ZernikeFirst+iZ]/_scaleFactor;
  }
  
  if (_printLevel>=2){
    std::cout << "DonutEngine: Parameters (Zernikes scaled) are = " << _nEle << " " << _rzero << " " << _bkgd << " " ;
    for (int iZ=0;iZ<nZernikeSize;iZ++){
      std::cout << _ZernikeArr[iZ] << " ";
    }
    std::cout << std::endl;
  }

  // first time through, initialize the _last_XXX variables
  // and set the counter to 0
  if (_first) {
    _last_nEle = 0.0;
    _last_rzero = 0.0;
    _last_bkgd = 0.0;
    _last_xDECam = 0.0;
    _last_yDECam = 0.0;
    _last_ZernikeArr.Dimension(nZernikeSize);
    _last_ZernikeArr.Activate();
    for (int iZ=0;iZ<nZernikeSize;iZ++){
      _last_ZernikeArr[iZ] = 0.0;
    }

    _stateCounter = 0;
    _first = false;
  } else {

    // reset State machine flags
    _statePupilMask = true;
    _statePupilFunc = true;
    _stateAtmos = true;

    // have any of the Pupil Wavefront parameters changed?
    for (int iZ=0;iZ<nZernikeSize;iZ++){
      if (_last_ZernikeArr[iZ] != _ZernikeArr[iZ]){
	_statePupilFunc = false;
      }
    }

    // check on state of usingDeltaWFM, if so then also set the state of the PupilFunc to false (needs recalculation)
    if (_usingDeltaWFM){
      _statePupilFunc = false;
    }


    // Pupil?
    if (_xDECam != _last_xDECam || _yDECam != _last_yDECam){
      _statePupilMask = false;
    }

    // Atmosphere?
    if (_last_rzero != _rzero){
      _stateAtmos = false;
    }

    if (_statePupilMask && _statePupilFunc && _stateAtmos){
      if (_debugFlag || _printLevel>=2 ){
	std::cout << "DonutEngine: no states have changed?" << std::endl;
      }
    }

    if (_printLevel>=2){
      std::cout << "DonutEngine: state machine flags " << _statePupilMask << " " << _statePupilFunc << " " << _stateAtmos << std::endl;
    }

  }

}

void DonutEngine::savePar(){
    
  // save parameter values for the next iteraiton
  _last_nEle = _nEle;  
  _last_rzero = _rzero;
  _last_bkgd = _bkgd;
  
  for (int iZ=0;iZ<nZernikeSize;iZ++){
    _last_ZernikeArr[iZ] = _ZernikeArr[iZ];
  }

  _last_xDECam = _xDECam;
  _last_yDECam = _yDECam;
  
}      

void DonutEngine::calcPupilMask(){

  clock_t start = clock();
    
  if (_printLevel>=2){
    std::cout << "DonutEngine: calcPupilMask x,y = " << _xDECam << " " << _yDECam << std::endl;
  }
  // input the PupilMask
  //     None:  generate the pupilMask from outer,inter Radius
  //     not implemented: string : readin the pupilMask from a fits file
  //     not implemented: numpy.ndarray: the numpy array with the pupilMask
  if (_inputPupilMask==""){

    // for DECam make pupil Mask with extra detail
    if (_iTelescope==0){
    
      // model all surfaces which contribute
      // only 3 obstructions matter: EndCapTop, Shroud, ChimneyAOut
      Real apertureRadius[3] = {0.301,0.298,0.264};
      Real apertureSlope[3] = {3.72e-4,2.31e-4,-1.79e-4};
      
      // calculate the offsets 
      // note that _xDECam,_yDECam is in [mm] not meters!!!
      Real offx[3];
      Real offy[3];
      for (int iap=0;iap<3;iap++){
	offx[iap] = _zemaxToDECamSignFlip * apertureSlope[iap] * _xDECam;
	offy[iap] = _zemaxToDECamSignFlip * apertureSlope[iap] * _yDECam;
      }

      //std::cout << "x,y,offx,offy = " << _xDECam << " " << _yDECam << " " << offx[0] << " " << offx[1] << " " << offx[2] << " " << offy[0] << " " << offy[1] << " " << offy[2] << " " << std::endl;

      // constants for the FilterExch
      Real filtExchSlope = 3.19e-4;
      Real filtExchBoxXmax = 0.167;
      Real filtExchBoxYmax = 0.320;
      Real filtExchArcP0 = 0.3685;
      Real filtExchArcP1 = -0.00559;
      Real filtExchArcP2 = -2.29;
      Real filtExchArcBoxXmax = 0.144;
      Real filtExchArcBoxYmax = 0.369;
      Real filtExchArcBoxYmin = 0.320;
      Real spiderWidth = 0.019 * (1462.526/4010.)/2.0;  // use 19mm thick
      //spiderWidth = 0.050 * (1462.526/4010.)/2.0;  // added May 15, 2014, try 50mm thick, actually works pretty well
      Real spiderSlope = 3.01e-4;

      // now loop over all bins, and build the pupil
      bool spiderMask(false);  //x1true mean pupil==1 (ie. clear)
      bool filtExchMask(false);
      bool annulusMask(false);
	
      for (int i=0;i<_nbin*_nbin;i++){

	// model the 3 circular apertures of the central obstruction
	double rhop[3] = {0.,0.,0.};
	for (int iap=0;iap<3;iap++){
	  double x = _xaxis(i) - offx[iap];
	  double y = _yaxis(i) - offy[iap];
	  rhop[iap] = sqrt(x*x+y*y)/apertureRadius[iap];	  	  
	}

	// pupil plane with central obscuration
	if ( _rho(i)<1.0 && rhop[0]>1.0 && rhop[1]>1.0 && rhop[2]>1.0 ){
	  annulusMask = true;
	} else {
	  annulusMask = false;
	}

	// model the Filter Exchange Mechanism
	double xx = fabs( _xaxis(i) -  _zemaxToDECamSignFlip * filtExchSlope * _xDECam);
	double yy = fabs( _yaxis(i) -  _zemaxToDECamSignFlip * filtExchSlope * _yDECam);
	
	if ( (xx<filtExchBoxXmax && yy<filtExchBoxYmax)  || 
	     (xx<filtExchArcBoxXmax && yy<filtExchArcBoxYmax && yy>=filtExchArcBoxYmin && yy < (filtExchArcP0+filtExchArcP1*xx+filtExchArcP2*xx*xx) ) ){	
	  filtExchMask = false;
	} else {
	  filtExchMask = true;
	}
    	
	// model the spider 
	// rotate by 45deg and make simple cut around spiderWidth
	double xxx = _xaxis(i) - _zemaxToDECamSignFlip * spiderSlope * _xDECam;
	double yyy = _yaxis(i) - _zemaxToDECamSignFlip * spiderSlope * _yDECam;
	double dxprime = xxx * cos(_M_PI_4) - yyy * sin(_M_PI_4);
	double dyprime = xxx * sin(_M_PI_4) + yyy * cos(_M_PI_4);
	
	if ( fabs(dxprime)<spiderWidth ||  fabs(dyprime)<spiderWidth ){
	  spiderMask = false;
	} else {
	  spiderMask = true;
	}
		
	// combine annulus, spider and filter Exchanger
	if (spiderMask && annulusMask && filtExchMask){
	  _pupilMask(i) = 1.0;
	} else {
	  _pupilMask(i) = 0.0;
	}
      } 

    }

  else if (_iTelescope==5){  //special case for DESI
    // build the pupil Mask
    Matrix dx(_nbin,_nbin);
    Matrix dy(_nbin,_nbin);
    dx = _xaxis;
    dy = _yaxis;
    
    Matrix rhoprime(_nbin,_nbin);
    Matrix dxprime(_nbin, _nbin);
    Matrix dyprime(_nbin, _nbin);
    Matrix _rhoellipse(_nbin, _nbin);

    Real spiderWidth = 0.0; //TBD
    Real invAsq = 0.57245381187;  //inverse square of the major axis of the aperture
    Real invBsq = 0.53424950221;  //inverse square of the minor axis of the aperture
    
    // create a grid for the aperture ellipse
    for (int i=0; i<_nbin * _nbin; i++){
      _rhoellipse(i) = sqrt(_xaxis(i)*_xaxis(i)*invAsq + _yaxis(i) * _yaxis(i) * invBsq) / _outerRadius;
    }
    
    Matrix spiderMask(_nbin,_nbin);
    Matrix annulusMask(_nbin,_nbin);

    Real a = 0.956343696  - 0.02034548701 * _xDECam + 0.05050383872 * _yDECam;
    Real b = 0.9664723832 - 0.0359776269  * _xDECam + 0.02464690834 * _yDECam;

    Real invasq = 1. / pow(a*1.37, 2);  // inverse square of obscuration major axis 
    Real invbsq = 1. / pow(b*1.37, 2);  // inverse square of obscuration minor axis 
    Real cy = -0.1233;   // center of obscuration ellipse on y axis
    Real cx = 0.13015;   // center of oscuration ellipse on x axis
    
    //calculate angle of rotation of obscuration ellipse
    Real phi = -353.750340246 + 192.46801126 * _xDECam + 112.10949469 * _yDECam;
    
    // calculate the slope and intercept of the vignetting feature, i.e. "bite"
    Real biteslope =  1.43295788    + 1.56       * _xDECam - 2.21165461 * _yDECam;
    Real biteint   = -8.81802987607 + 2.71588835 * _xDECam + 3.84451752 * _yDECam;

    // Declare a few variables we need for the mask loop below
    Real lhs;
    Real rhs;
    //Real bitex;
    Real bitey;

    for (int i=0;i<_nbin*_nbin;i++){
        // build a mask for the obscuration ellipse
	lhs = cos(phi) * (dx(i) - cx) + (dy(i) - cy) * sin(phi);
	rhs = sin(phi) * (dx(i) - cx) - (dy(i) - cy) * cos(phi);
	rhoprime(i) = sqrt(lhs * lhs * invasq + rhs * rhs * invbsq) / _innerRadius;
	
	// rotate by 45deg, cut around spiderWidth, rotate back
	dxprime(i) = dx(i) * cos(_M_PI_4) - dy(i) * sin(_M_PI_4);
	dyprime(i) = dx(i) * sin(_M_PI_4) + dy(i) * cos(_M_PI_4);
	

	if ( fabs(dxprime(i))<spiderWidth ||  fabs(dyprime(i))<spiderWidth ){
	  spiderMask(i) = 0.0;
	} else {
	  spiderMask(i) = 1.0;
	}
	
	// pupil plane with central obscuration
	if ( _rhoellipse(i)<1.0 && rhoprime(i)>1. ){
	  annulusMask(i) = 1.0;
	} else {
	  annulusMask(i) = 0.0;
	}
	// calculate the bite function and incorporate it into the pupil mask
        bitey = biteint + biteslope * dx(i);
	if (dy(i) < bitey){
	  annulusMask(i) = 0.0;
	}
	
	// combine spider and annulus
	if (spiderMask(i)==1.0 && annulusMask(i)==1.0){
	  _pupilMask(i) = 1.0;
	} else {
	  _pupilMask(i) = 0.0;
	}
    } 
    
    
  }
  else {
      // input parameters will define pupilMask
      // parameters: spiderWidth, innerRadius
      
      // build the pupil Mask
      Matrix dx(_nbin,_nbin);
      Matrix dy(_nbin,_nbin);
      dx = _xaxis;
      dy = _yaxis;
      
      Matrix rhoprime(_nbin,_nbin);
      Matrix dxprime(_nbin,_nbin);
      Matrix dyprime(_nbin,_nbin);
      //Real spiderWidth = 0.0028/2.0; 
      Real spiderWidth = 0.0;   //Kludge on Sept 22, 2014
      
      Matrix spiderMask(_nbin,_nbin);
      Matrix annulusMask(_nbin,_nbin);
      
      for (int i=0;i<_nbin*_nbin;i++){
	rhoprime(i) = sqrt(dx(i)*dx(i) + dy(i)*dy(i))/_innerRadius;
	
	// rotate by 45deg, cut around spiderWidth, rotate back
	dxprime(i) = dx(i) * cos(_M_PI_4) - dy(i) * sin(_M_PI_4);
	dyprime(i) = dx(i) * sin(_M_PI_4) + dy(i) * cos(_M_PI_4);
	
	if ( fabs(dxprime(i))<spiderWidth ||  fabs(dyprime(i))<spiderWidth ){
	  spiderMask(i) = 0.0;
	} else {
	  spiderMask(i) = 1.0;
	}
	
	// pupil plane with central obscuration
	if ( _rho(i)<1.0 && rhoprime(i)>1. ){
	  annulusMask(i) = 1.0;
	} else {
	  annulusMask(i) = 0.0;
	}
	
	// combine spider and annulus
	if (spiderMask(i)==1.0 && annulusMask(i)==1.0){
	  _pupilMask(i) = 1.0;
	} else {
	  _pupilMask(i) = 0.0;
	}
      } 

    } // test for iTelescope==0

  } else {
    // inputPupilMask - array defining the mask, instead of using parameters
    // DISABLED NOW!!!
    std::cout << "DonutEngine:  ERROR inputPupilMask is currently disabled " << std::endl;
  }

  // save sqrt of the pupilMask.sum()
  _pupilSNorm = 0.0;
  for (int i=0;i<_nbin*_nbin;i++){
    _pupilSNorm += _pupilMask(i);
  }    
  _pupilSNorm = sqrt(_pupilSNorm);
  
  clock_t stop = clock();
  _timePupilMask += (stop-start)/(Real)CLOCKS_PER_SEC;

}

void DonutEngine::calcPupilFunc(){

  clock_t start = clock();

  if (_printLevel>=2){
    std::cout << "DonutEngine: calcPupilFunc " << std::endl;
  }
  
  // Zernike terms
  for (int iZ=0;iZ<nZernikeSize;iZ++){
    if (_ZernikeArr[iZ] != _last_ZernikeArr[iZ]){

      Matrix zernikeTemp =  _zernikeObject->_zernikeTerm[iZ+1]; // need [iZ+1] since ZernikeTerm includes the Piston term
      for (int i=0;i<_nbin*_nbin;i++){
	_pupilWaveZernike(i) = _pupilWaveZernike(i) +  ( _ZernikeArr[iZ] - _last_ZernikeArr[iZ] ) * zernikeTemp(i);
      }
      
    }
  }

  //add deltaWFM if desired
  if (_usingDeltaWFM){
     for (int i=0;i<_nbin*_nbin;i++){
       _pupilWaveZernikePlusDelta(i) = _pupilWaveZernike(i) + _deltaWFM(i);
     }
  }
    
  // calculate the pupilFunc(tion) from the pupilWave(front)
  //Complex I = Complex(0.0,1.0);
  Complex twopiI = Complex(0.0,2.0*_M_PI);
  for (int i=0;i<_nbin*_nbin;i++){
    if (_pupilMask(i)==0.0){
      _pupilFunc(i) = 0.0;
      _pupilFuncStar(i) = 0.0;
    } else {
      if (_usingDeltaWFM){
	_pupilFunc(i) = _pupilMask(i) * exp(twopiI  *_pupilWaveZernikePlusDelta(i));   // no lambda here, so units are in waveLength
      } else{
	_pupilFunc(i) = _pupilMask(i) * exp(twopiI  *_pupilWaveZernike(i));   // no lambda here, so units are in waveLength
      }
      _pupilFuncStar(i) = conj(_pupilFunc(i));
    }
  }

  if (_debugFlag && nCallsCalcAll<=1){    
    Matrix pupil(_nbin,_nbin);
    for (int i=0;i<_nbin*_nbin;i++){
      pupil(i) = _pupilWaveZernike(i) * _pupilMask(i);
    }
    toFits(_fptr,pupil);
    toFits(_fptr,_pupilFunc);
  }

  clock_t stop = clock();
  _timePupilFunc += (stop-start)/(Real)CLOCKS_PER_SEC;

}
 
        
void DonutEngine::calcOptics(){
           
  clock_t start = clock();

  if (_printLevel>=2){
    std::cout << "DonutEngine: calcOptics " << std::endl;
  }
                
  // calculate the PSF from the Pupil function (use ifft to match Zemax output! )
  MatrixC ishftpupilFunc(_nbin,_nbin);
  ishftpupilFunc = _pupilFunc;
  fftShift(ishftpupilFunc); //was InvShift, but these are the same as long as _nbin is even
  _ifftInputArray = ishftpupilFunc;  
  _ifft2PlanC->execute();

  // don't shift G and G*, only psfOptics   (normalize to sqrt(Area*NbinsTotal))
  Real normalizationG = 1.0/(_nbin*_pupilSNorm);
  Matrix unshftpsfOptics(_nbin,_nbin);
  for (int i=0;i<_nbin*_nbin;i++){
    _calcG(i) = _ifftOutputArray(i) * normalizationG;
    _calcGstar(i) = conj(_calcG(i));
    unshftpsfOptics(i) = real(_calcG(i)*_calcGstar(i));
  }

  if (_debugFlag  && nCallsCalcAll<=1){
    _psfOptics = unshftpsfOptics;
    fftShift(_psfOptics);

    toFits(_fptr,_psfOptics);
  }

  // now take the Fourier Transform of the PSF for use in Convolution
  // Q: is an fft of an unshifted fft shifted?
  realToComplex(unshftpsfOptics,_fftInputArray);
  _fft2PlanC->execute();
  _ftsOptics = _fftOutputArray;  

  //_fftrtcInputArray = unshftpsfOptics;
  //_fft2rtcPlanC->execute();
  //_ftsOptics = _fftrtcOutputArray;  

  // for (int k=0;k<512*512;k++){
  //   if (abs(_ftsOptics(k)-_fftrtcOutputArray(k))>1e-10){
  //     std::cout << k << " " << _ftsOptics(k) << " " << _fftrtcOutputArray(k) << std::endl;
  //   }
  // }

  clock_t stop = clock();
  _timeOptics += (stop-start)/(Real)CLOCKS_PER_SEC;

}

void DonutEngine::shiftnormG(){
  // shift and normalize Magnitude of _calcG, also get its Phase
  fftShift(_calcG);
  Real summagG(0.);
  for (int i=0;i<_nbin*_nbin;i++){
    _magG(i) = abs(_calcG(i));
    _phaseG(i) = arg(_calcG(i));
    summagG = summagG + _magG(i);
  }
  for (int i=0;i<_nbin*_nbin;i++){
    _magG(i) = _magG(i)/summagG;
  }
}
    

void DonutEngine::calcAtmos(){

  clock_t start = clock();

  if (_printLevel>=2){
    std::cout << "DonutEngine: calcAtmos" << std::endl;
  }
  
  // calculate Kolmogorov dist, use shifted radius Array, instead of shifting this everytime!
  Matrix shftarrAtmos(_nbin,_nbin);
  Real fivethirds(5./3.);
  Real shftarrAtmosMax(0.);
  for (int i=0;i<_nbin*_nbin;i++){  
    shftarrAtmos(i) = exp(-3.44*pow(_shftrAtmos(i)*_waveLength*_fLength/_rzero,fivethirds));
    if (shftarrAtmos(i)>shftarrAtmosMax){
      shftarrAtmosMax = shftarrAtmos(i);
    }
  }
  // normalize (for now) to match python code
  shftarrAtmos /= shftarrAtmosMax;

  // take the inverse FT to get the Atmosphere's PSF
  realToComplex(shftarrAtmos,_ifftInputArray);  
  _ifft2PlanC->execute();

  //_fftrtcInputArray = shftarrAtmos;
  //_fft2rtcPlanC->execute();

  Matrix unshftpsfAtmos(_nbin,_nbin);
  // normalization of unshftpsfAtmos is very close to the maximum value divided by _nbin*_nbin
  // but is a few percent off from that - so just normalize so the sum==1.0 
  Real atmosNormalization(0.);
  for (int i=0;i<_nbin*_nbin;i++){  
    unshftpsfAtmos(i) = abs(_ifftOutputArray(i))/(_nbin*_nbin);
    //unshftpsfAtmos(i) = abs(_fftrtcOutputArray(i))/(_nbin*_nbin);
    atmosNormalization += unshftpsfAtmos(i);
  }
  unshftpsfAtmos *= (1.0/atmosNormalization);

  fftShift(unshftpsfAtmos,_psfAtmos);  //do by default now...

  if (_debugFlag  && nCallsCalcAll<=1){
    fftShift(unshftpsfAtmos,_psfAtmos);
    toFits(_fptr,_psfAtmos);
  }

                                
  // to convolve with the other sources of PSF, take the FT of the Atmosphere's PSF
  realToComplex(unshftpsfAtmos,_fftInputArray); 
  _fft2PlanC->execute();
  _ftsAtmos = _fftOutputArray; 

  //_fftrtcInputArray = unshftpsfAtmos;
  //_fft2rtcPlanC->execute();
  //_ftsAtmos = _fftrtcOutputArray; 

  clock_t stop = clock();
  _timeAtmos += (stop-start)/(Real)CLOCKS_PER_SEC;

}

void DonutEngine::calcFTPixel(){

  // define the pixelBox
  _pixelBox.Dimension(_nbin,_nbin);
  _pixelBox.Activate();
  Real sumOfBox(0.);
  for (int i=0;i<_nbin*_nbin;i++){
    if ( fabs(_xpsf(i)) <= _pixelSize/2.0 &&  fabs(_ypsf(i)) <= _pixelSize/2.0 ){
      _pixelBox(i) = 1.0;
      sumOfBox += abs(_pixelBox(i));
    } else {
      _pixelBox(i) = 0.0;
    }
  }

  if (_debugFlag  && nCallsCalcAll<=1){
    toFits(_fptr,_pixelBox);
  }

  // check if grid is too sparse, then just set the FT(pixel) to all ones
  if (sumOfBox>0.){
    _pixelBox *= (1.0/sumOfBox);           // normalize to 1.0
    fftShift(_pixelBox,_fftInputArray);  // was InvShift, also don't need pixelBox anymore - this shifts in place
    _fft2PlanC->execute();
    _ftsPixel = _fftOutputArray; 
  } else {
    _ftsPixel = 1.0;  // sets whole array to 1.0
  }
            

  if (_debugFlag  && nCallsCalcAll<=1){
    toFits(_fptr,_ftsPixel);
  }

}

void DonutEngine::calcConvolute(){

  clock_t start = clock();

  if (_printLevel>=2){
    std::cout << "DonutEngine: calcConv" << std::endl;
  }
        
  // convolution (now doing it as F-1{F(Optics) F(Atmos) F(Pixels)
  MatrixC prodFts(_nbin,_nbin);
  for (int i=0;i<_nbin*_nbin;i++){
    prodFts(i) = _ftsOptics(i) * _ftsAtmos(i) * _ftsPixel(i);
  }

  _ifftInputArray = prodFts;
  _ifft2PlanC->execute();
  fftShift(_ifftOutputArray);

  // normalize and take absolute value
  Real nsqNorm = 1.0/(_nbin*_nbin);

  // save calculated image
  for (int i=0;i<_nbin*_nbin;i++){
    _convOpticsAtmosPixel(i) = abs(_ifftOutputArray(i)) * nsqNorm;
  }

  if (_debugFlag  && nCallsCalcAll<=1){
    toFits(_fptr,_convOpticsAtmosPixel);
  }

  clock_t stop = clock();
  _timeConvolute += (stop-start)/(Real)CLOCKS_PER_SEC;

}
        
            
void DonutEngine::calcPixelate(){

  clock_t start = clock();

  if (_printLevel>=2){
    std::cout << "DonutEngine: calcPixelate" << std::endl;
  }
  
  // integer value of ngridperPixel?
  bool useNgridperPixel = (int(_ngridperPixel) == _ngridperPixel);

  if (useNgridperPixel){
    // pixelate - now by just selecting every nth grid point!
    // Q: is this there where we are getting offset from the grid center?  should we start at the stride/2 grid point and
    //    then go every stride grid point???
    // A: actually, this may be where the grid is offset from 0, but there is no solution with an even number of grid points
    //    and 0.0,0.0 in between grid points and ngridperPixel an even value,  this could work if ngridperPixel is odd!!
    int index(0);
    int stride = (int) _ngridperPixel;
    int pixIndex(0);
    Real gridNorm = _ngridperPixel*_ngridperPixel;
    for (int j=0;j<_nbin;j=j+stride){
      for (int i=0;i<_nbin;i=i+stride){

	_valPixelCenters(pixIndex) = _convOpticsAtmosPixel(index) * gridNorm;


	//cout << i << " " << j << " " <<  index << " " << pixIndex << " " << _convOpticsAtmosPixel(index) << " " << _valPixelCenters(pixIndex) << endl;
	index = index + stride;
	pixIndex++;

	// numpy equivalent
	//  _xPixels = _xpsf[0:_nbin:_ngridperPixel,0:_nbin:_ngridperPixel]
	//  _yPixels = _ypsf[0:_nbin:_ngridperPixel,0:_nbin:_ngridperPixel]

      }
      index = index + _nbin*(stride-1);
    }

  } else{  // use a GSL routine
    // pixelate with interpolation
    //_xPixelCenters,_yPixelCenters,_valPixelCenters = _getPixelated(0.0,0.0,_pixelSize,_nPixels,_xpsf,_ypsf,_convOpticsAtmosPixel);
    // normalize
    //_valPixelCenters = _valPixelCenters/_valPixelCenters.sum();
    
  }

  if (_debugFlag  && nCallsCalcAll<=1){
    toFits(_fptr,_valPixelCenters);
  }

  clock_t stop = clock();
  _timePixelate += (stop-start)/(Real)CLOCKS_PER_SEC;

}
 

    // void getPixelated(self,xPixOffset,yPixOffset,pixelSize,nPixels,xpsf,ypsf,convPSF):

    //     // xPixOffset and yPixOffset go from -1 to 1 inside a pixel
    //     xLocOffset = xPixOffset*pixelSize/2.0
    //     yLocOffset = yPixOffset*pixelSize/2.0   // from pixel to local coordinates
        
    //     // setup interpolation, and get values at pixel centers given that the central pixel is at a certain position
    //     // get nPixels x nPixels sized array
        
    //     xLo = -(nPixels/2)*pixelSize + xLocOffset
    //     xHi = (nPixels/2)*pixelSize + xLocOffset
    //     yLo = -(nPixels/2)*pixelSize + yLocOffset
    //     yHi = (nPixels/2)*pixelSize + yLocOffset
    //     yPixelCenters,xPixelCenters = itricks.mgrid[yLo:yHi:(nPixels)*1j,xLo:xHi:(nPixels)*1j]
        
    //     myInterp = interpolate.RectBivariateSpline(ypsf[:,0],xpsf[0,:],convPSF,kx=3,ky=3)
    //     valPixelCenters = myInterp(yPixelCenters[:,0],xPixelCenters[0,:])
        
    //     return xPixelCenters,yPixelCenters,valPixelCenters


void DonutEngine::calcDerivatives(double* image, int ny, int nx, 
				  double* weight, int my, int mx){
  calcDerivatives(image,weight);
}

void DonutEngine::calcDerivatives(Real* image, Real* weight){

  nCallsCalcDerivative++;

  clock_t start = clock();

  // calculate derivatives, put in the _dChi2dpar array

  // needed arrays:  calcImage, ftsAtmos, ftsPixel, calcG, calcGstar
  // local arrays (save time by avoiding re-reserving memory?):
  //     Qpixels, Q, Qtilde, QQ, QQtilde, QQQ QQQtilde

  // arrays we will need
  Matrix Qpixels(_nPixels,_nPixels);
  Matrix Q(_nbin,_nbin);
  Q = 0.0;
  MatrixC Qtilde(_nbin,_nbin);
  MatrixC QQ(_nbin,_nbin);
  MatrixC QQtilde(_nbin,_nbin);
  MatrixC QQQ(_nbin,_nbin);
  MatrixC QQQtilde(_nbin,_nbin);

  clock_t stop = clock();
  _timeDerivatives0 += (stop-start)/(Real)CLOCKS_PER_SEC;

  start = clock();

  // calculate Q = W(I-N)
  for (int i=0;i<_nPixels*_nPixels;i++){
    Qpixels(i) = weight[i] * (_calcImage(i) - image[i]);
  }

  // pepper the image on a full-size grid
  int index(0);
  int pixIndex(0);
  int stride = (int) _ngridperPixel;

  for (int j=0;j<_nbin;j=j+stride){
    for (int i=0;i<_nbin;i=i+stride){
	Q(index) = Qpixels(pixIndex);

	index = index + stride;
	pixIndex++;
    }
    index = index + _nbin*(stride-1);
  }

  // FT^-1{Q}
  fftShift(Q); //was InvShift
  realToComplex(Q,_ifftInputArray);  
  _ifft2PlanC->execute();
  Qtilde = _ifftOutputArray;

  //_fftrtcInputArray = Q;
  //_fft2rtcPlanC->execute();
  //Qtilde = _fftrtcOutputArray;


  // F{ Qtilde * ftsAtmos * ftsPixel}
  for (int i=0;i<_nbin*_nbin;i++){
    QQ(i) = Qtilde(i) * _ftsAtmos(i) * _ftsPixel(i);
  }
  _fftInputArray = QQ;
  _fft2PlanC->execute();
  Real QQnorm = 1.0/(_nbin * _nbin);
  QQtilde = _fftOutputArray;
  //QQtilde *= QQnorm;

  // QQQstartilde and QQQtilde are complex conj, can use that, instead of separate calcs
  // since QQtilde is all real
  // F{ G * QQtilde } 

  for (int i=0;i<_nbin*_nbin;i++){ 
    QQQ(i) = _calcG(i) * QQtilde(i);
  }
  _fftInputArray = QQQ;
  _fft2PlanC->execute();
  Real QQQnorm = 1.0/(_nbin * _nbin);
  fftShift(_fftOutputArray);
  QQQtilde = _fftOutputArray;
  Real QQandQQQnorm = QQnorm * QQQnorm;
  QQQtilde *= QQandQQQnorm;

  stop = clock();
  _timeDerivatives1 += (stop-start)/(Real)CLOCKS_PER_SEC;

  start = clock();

  // dg*(x)/dalpha * QQQtilde
  Vector dChi2dzern(nZernikeSize);
  
  // Loop over Zernike terms - only need nZernikeSize=nZernikeTerm-1 of them!!!
  Complex minustwopiI(0.0,-2.0*_M_PI);
  Real minustwopi(-2.0*_M_PI);
  for (int iZ=0;iZ<nZernikeSize;iZ++){

    MatrixC dgdalphaStar(_nbin,_nbin);
    dgdalphaStar = 0.0;
    Matrix zernikeTemp = _zernikeObject->_zernikeTerm[iZ+1];  // since nZernikeSize = nZernikeTerm-1
    for (int i=0;i<_nbin*_nbin;i++){
      // these are complex conj too, can calculate more compactly, 4.0 * Re{} below
      // move 2piI below, was
      // dgdalphaStar(i) = minustwopiI *  _pupilFuncStar(i) * zernikeTemp(i);
      dgdalphaStar(i) =  _pupilFuncStar(i) * zernikeTemp(i);
      ////dgdalpha = 2.0 * numpy.pi * 1j *  _pupilFunc * _zernikeObject->_zernikeTerm[iZ]
    }

    ////dChi2dalpha[iZ] = 2.0 * _nEle *  (dgdalphaStar * QQQtilde).sum()  + 2.0 * _nEle *  (dgdalpha * QQQstartilde).sum()
    ////is equal to    dChi2dzern[iZ] = 4.0 * _nEle *  ((dgdalphaStar * QQQtilde).real).sum()
    dChi2dzern[iZ] = 0.0;
    for (int i=0;i<_nbin*_nbin;i++){
      //      dChi2dzern[iZ] += real(dgdalphaStar(i)*QQQtilde(i));
      dChi2dzern[iZ] -= imag(dgdalphaStar(i)*QQQtilde(i));  //note its -= not += !!
    }
    dChi2dzern[iZ] *= (4.0 * _nEle * minustwopi * 86.8692) ; 
    // why oh why am I off by this weird number??!!
    //    dChi2dzern[iZ] *= 86.8692;

    // need extra scaling  with scaleFactor - 3 factors, 1 for Zernikes, 2 for grid
    dChi2dzern[iZ] = dChi2dzern[iZ]/(_scaleFactor*_scaleFactor*_scaleFactor);

  }

  // also calculate derivative of Nele,bckg and rzero
  Real dChi2dnele(0.);
  Real dChi2dbkgd(0.);
  for (int i=0;i<_nPixels*_nPixels;i++){
    dChi2dnele += 2.0 * Qpixels(i) * (_calcImage(i) - _bkgd) / _nEle; 
    dChi2dbkgd += 2.0 * Qpixels(i);
  }

  // calculate the derivative of rzero, if desired
  // assumes that calcAll has been called prior to calling calcDerivative!!!
  Real dChi2drzero(0.);
  if (_calcRzeroDerivative){

    // calculate Chi2 value now
    Real chi2nominal(0.);
    for (int i=0;i<_nPixels*_nPixels;i++){
      chi2nominal = chi2nominal + weight[i] * (_calcImage(i) - image[i]) * (_calcImage(i) - image[i]);
    }

    // set the x,y position for calculating the pupil
    _anotherDonutEngine->setXYDECam(_xDECam,_yDECam);

    // calculate Chi2 value at rzero + delta
    Real delta(0.0005);

    // calculate chi2 at +delta
    Vector parToUse;
    parToUse.Dimension(npar);
    parToUse.Activate();
    // save input parameters
    for (int ipar=0;ipar<npar;ipar++){
      parToUse[ipar] = _parCurrent[ipar];
    }
    parToUse[ipar_rzero] = parToUse[ipar_rzero] + delta;

    _anotherDonutEngine->calcAll(parToUse);

    // calculate Chi2 value at this +delta value
    Real chi2deltaplus(0.);
    for (int i=0;i<_nPixels*_nPixels;i++){
      chi2deltaplus = chi2deltaplus + weight[i] * pow(_anotherDonutEngine->getImage()(i) - image[i],2);
    }

    // calculate chi2 at -delta
    parToUse[ipar_rzero] = parToUse[ipar_rzero] - 2. * delta;

    _anotherDonutEngine->calcAll(parToUse);

    // calculate Chi2 value at this +delta value
    Real chi2deltaminus(0.);
    for (int i=0;i<_nPixels*_nPixels;i++){
      chi2deltaminus = chi2deltaminus + weight[i] * pow(_anotherDonutEngine->getImage()(i) - image[i],2);
    }

    dChi2drzero = (chi2deltaplus - chi2deltaminus)/(2.*delta);
  }


  // fill return array
  _dChi2dpar[ipar_nEle] = dChi2dnele;
  _dChi2dpar[ipar_bkgd] = dChi2dbkgd;
  _dChi2dpar[ipar_rzero] = dChi2drzero;
  for (int iZ=0;iZ<nZernikeSize;iZ++){
    _dChi2dpar[ipar_ZernikeFirst+iZ] = dChi2dzern[iZ];
  }

  // print out Derivatives
  if (_printLevel>=2){
    std::cout << "DonutEngine: Derivatives are = " ;
    for (int ipar=0;ipar<npar;ipar++){
      std::cout << _dChi2dpar[ipar] << " ";
    }
    std::cout << std::endl;
  }


  stop = clock();
  _timeDerivatives2 += (stop-start)/(Real)CLOCKS_PER_SEC;

}


    // void calcImageDerivatives(self,paramArray,paramErrArray,fixParamArray,image,weight):
    //     // calc dI(u)/dalpha for all parameters (even r0) and return maxtrix

    //     // calculate starting image
    //     _calcAll(paramArray)
    //     calcImageSaved = _calcImage
        
    //     // transfer matrix A[iPixel,iPar]
    //     transArray = matrix(numpy.zeros((_nPixels*_nPixels,npar)))

    //     // save current parameter values
    //     parSaved = _parCurrent.copy()
    //     calcMatrix = matrix(_calcImage)
    //     calcMatrix.resize(_nPixels*_nPixels)

    //     // loop over desired parameters, change them one at a time and get their derivative per pixel
    //     for ipar in range(npar):
    //         if fixParamArray[ipar]==0:
    //             parDeriv = copy.deepcopy(parSaved)
    //             delta = paramErrArray[ipar]
    //             if delta<=0.0:
    //                 delta = 1.0e-4
    //             parDeriv[ipar] = parDeriv[ipar] + delta
    //             _calcAll(parDeriv)
    //             thisMatrix = matrix(_calcImage)
    //             thisMatrix.resize(_nPixels*_nPixels)                
    //             deltaImage = (thisMatrix - calcMatrix)/delta
    //             transArray[:,ipar] = deltaImage.transpose()

    //     // restore original values and image
    //     _parCurrent = parSaved.copy()
    //     _calcAll(_parCurrent)

    //     // return the matrix
    //     return transArray,calcImageSaved
    

void DonutEngine::realToComplex(Matrix& in, MatrixC& out){

  int n = in.Nx() * in.Ny();
  // assume both arrays have same size 

  for (int i=0;i<n;i++){
    out(i) = (Complex) in(i);
  }

}


void DonutEngine::absComplex(MatrixC& in, Matrix& out){

  int nx = in.Nx();
  int ny = in.Nx();

  for (int i=0;i<nx*ny;i++){
    out(i) = abs(in(i));
  }

}


void DonutEngine::resetTimers(){
  _timePupilMask=0.0;
  _timePupilFunc= 0.0;
  _timeOptics= 0.0;
  _timeAtmos= 0.0;
  _timeConvolute= 0.0;
  _timePixelate= 0.0;
  _timeDerivatives0 = 0.0;
  _timeDerivatives1 = 0.0;
  _timeDerivatives2 = 0.0;
}

void DonutEngine::printTimers(){

  std::cout << "DonutEngine Timers " << std::endl;
  std::cout << "     Pupil Mask     = " << _timePupilMask << std::endl;
  std::cout << "     Pupil Func     = " << _timePupilFunc << std::endl;
  std::cout << "     Optics         = " << _timeOptics << std::endl;
  std::cout << "     Atmos          = " << _timeAtmos << std::endl;
  std::cout << "     Convolute      = " << _timeConvolute << std::endl;
  std::cout << "     Pixelate       = " << _timePixelate << std::endl;
  std::cout << "     Derivatives0   = " << _timeDerivatives0 << std::endl;
  std::cout << "     Derivatives1   = " << _timeDerivatives1 << std::endl;
  std::cout << "     Derivatives2   = " << _timeDerivatives2 << std::endl;
}

void DonutEngine::getvXaxis(double** ARGOUTVIEW_ARRAY2, int* DIM1, int* DIM2){
  *DIM1 = _xaxis.Ny();
  *DIM2 = _xaxis.Nx();
  *ARGOUTVIEW_ARRAY2 = _xaxis();
}

void DonutEngine::getvYaxis(double** ARGOUTVIEW_ARRAY2, int* DIM1, int* DIM2){
  *DIM1 = _yaxis.Ny();
  *DIM2 = _yaxis.Nx();
  *ARGOUTVIEW_ARRAY2 = _yaxis();
}

void DonutEngine::getvRho(double** ARGOUTVIEW_ARRAY2, int* DIM1, int* DIM2){
  *DIM1 = _rho.Ny();
  *DIM2 = _rho.Nx();
  *ARGOUTVIEW_ARRAY2 = _rho();
}

void DonutEngine::getvTheta(double** ARGOUTVIEW_ARRAY2, int* DIM1, int* DIM2){
  *DIM1 = _theta.Ny();
  *DIM2 = _theta.Nx();
  *ARGOUTVIEW_ARRAY2 = _theta();
}

void DonutEngine::getvPixelBox(Complex** ARGOUTVIEW_ARRAY2, int* DIM1, int* DIM2){
  *DIM1 = _pixelBox.Ny();
  *DIM2 = _pixelBox.Nx();
  *ARGOUTVIEW_ARRAY2 = _pixelBox();
}

void DonutEngine::getvFtsPixel(Complex** ARGOUTVIEW_ARRAY2, int* DIM1, int* DIM2){
  *DIM1 = _ftsPixel.Ny();
  *DIM2 = _ftsPixel.Nx();
  *ARGOUTVIEW_ARRAY2 = _ftsPixel();
}

void DonutEngine::getvPupilWaveZernike(double** ARGOUTVIEW_ARRAY2, int* DIM1, int* DIM2){
  *DIM1 = _pupilWaveZernike.Ny();
  *DIM2 = _pupilWaveZernike.Nx();
  *ARGOUTVIEW_ARRAY2 = _pupilWaveZernike();
}

void DonutEngine::getvPupilWaveZernikePlusDelta(double** ARGOUTVIEW_ARRAY2, int* DIM1, int* DIM2){
  *DIM1 = _pupilWaveZernike.Ny();
  *DIM2 = _pupilWaveZernike.Nx();
  *ARGOUTVIEW_ARRAY2 = _pupilWaveZernikePlusDelta();
}

void DonutEngine::getvDeltaWFM(double** ARGOUTVIEW_ARRAY2, int* DIM1, int* DIM2){
  *DIM1 = _deltaWFM.Ny();
  *DIM2 = _deltaWFM.Nx();
  *ARGOUTVIEW_ARRAY2 = _deltaWFM();
}

void DonutEngine::getvPupilMask(double** ARGOUTVIEW_ARRAY2, int* DIM1, int* DIM2){
  *DIM1 = _pupilMask.Ny();
  *DIM2 = _pupilMask.Nx();
  *ARGOUTVIEW_ARRAY2 = _pupilMask();
}

void DonutEngine::getvPupilFunc(Complex** ARGOUTVIEW_ARRAY2, int* DIM1, int* DIM2){
  *DIM1 = _pupilFunc.Ny();
  *DIM2 = _pupilFunc.Nx();
  *ARGOUTVIEW_ARRAY2 = _pupilFunc();
}

void DonutEngine::getvMagG(double** ARGOUTVIEW_ARRAY2, int* DIM1, int* DIM2){
  *DIM1 = _magG.Ny();
  *DIM2 = _magG.Nx();
  *ARGOUTVIEW_ARRAY2 = _magG();
}

void DonutEngine::getvPhaseG(double** ARGOUTVIEW_ARRAY2, int* DIM1, int* DIM2){
  *DIM1 = _phaseG.Ny();
  *DIM2 = _phaseG.Nx();
  *ARGOUTVIEW_ARRAY2 = _phaseG();
}

void DonutEngine::getvPsfOptics(double** ARGOUTVIEW_ARRAY2, int* DIM1, int* DIM2){
  *DIM1 = _psfOptics.Ny();
  *DIM2 = _psfOptics.Nx();
  *ARGOUTVIEW_ARRAY2 = _psfOptics();
}

void DonutEngine::getvFtsOptics(Complex** ARGOUTVIEW_ARRAY2, int* DIM1, int* DIM2){
  *DIM1 = _ftsOptics.Ny();
  *DIM2 = _ftsOptics.Nx();
  *ARGOUTVIEW_ARRAY2 = _ftsOptics();
}

void DonutEngine::getvPsfAtmos(double** ARGOUTVIEW_ARRAY2, int* DIM1, int* DIM2){
  *DIM1 = _psfAtmos.Ny();
  *DIM2 = _psfAtmos.Nx();
  *ARGOUTVIEW_ARRAY2 = _psfAtmos();
}

void DonutEngine::getvFtsAtmos(Complex** ARGOUTVIEW_ARRAY2, int* DIM1, int* DIM2){
  *DIM1 = _ftsAtmos.Ny();
  *DIM2 = _ftsAtmos.Nx();
  *ARGOUTVIEW_ARRAY2 = _ftsAtmos();
}

void DonutEngine::getvConvOpticsAtmosPixel(double** ARGOUTVIEW_ARRAY2, int* DIM1, int* DIM2){
  *DIM1 = _convOpticsAtmosPixel.Ny();
  *DIM2 = _convOpticsAtmosPixel.Nx();
  *ARGOUTVIEW_ARRAY2 = _convOpticsAtmosPixel();
}

void DonutEngine::getvValPixelCenters(double** ARGOUTVIEW_ARRAY2, int* DIM1, int* DIM2){
  *DIM1 = _valPixelCenters.Ny();
  *DIM2 = _valPixelCenters.Nx();
  *ARGOUTVIEW_ARRAY2 = _valPixelCenters();
}

void DonutEngine::getvImage(double** ARGOUTVIEW_ARRAY2, int* DIM1, int* DIM2){
  *DIM1 = _calcImage.Ny();
  *DIM2 = _calcImage.Nx();
  *ARGOUTVIEW_ARRAY2 = _calcImage();
}

void DonutEngine::getParCurrent(double** ARGOUTVIEW_ARRAY1, int* DIM1){
  *DIM1 = _parCurrent.Nx();
  *ARGOUTVIEW_ARRAY1 = _parCurrent;
}

void DonutEngine::getDerivatives(double** ARGOUTVIEW_ARRAY1, int* DIM1){
  *DIM1 = _dChi2dpar.Nx();
  *ARGOUTVIEW_ARRAY1 = _dChi2dpar;
}


void DonutEngine::toFits(fitsfile* fptr, MatrixC& in){

  /* Create the image */
  long naxes[2];
  naxes[0] = in.Nx();
  naxes[1] = in.Ny();
  long Nsq = naxes[0] * naxes[1];
  long fpixel=1;
  int status = 0;         /* initialize status before calling fitsio routines */

  if (fits_create_img(fptr,DOUBLE_IMG,2,naxes,&status)){
    fits_report_error(stderr, status);
  }

  /* Write the array of doubles to the image */
  Matrix temp(naxes[1],naxes[0]);
  for (int i=0;i<Nsq;i++){
    temp(i) = abs(in(i));
  }
  if (fits_write_img(fptr, TDOUBLE, fpixel, Nsq, temp(), &status)){
    fits_report_error(stderr, status);
  }

}


void DonutEngine::toFits(fitsfile* fptr, Matrix& in){

  /* Create the image */
  long naxes[2];
  naxes[0] = in.Nx();
  naxes[1] = in.Ny();
  long Nsq = naxes[0] * naxes[1];
  long fpixel=1;
  int status = 0;         /* initialize status before calling fitsio routines */

  if (fits_create_img(fptr,DOUBLE_IMG,2,naxes,&status)){
    fits_report_error(stderr, status);
  }

  /* Write the array of doubles to the image */
  Real* ptrIn = in();
  if (fits_write_img(fptr, TDOUBLE, fpixel, Nsq, ptrIn, &status)){
    fits_report_error(stderr, status);
  }

}




