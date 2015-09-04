//
// DonutEngine.h:  Engine to calculate focal plane image from 
//                  pupil plane Zernike expansion
// Copyright (C) 2011 Aaron J. Roodman, SLAC National Accelerator Laboratory, Stanford University
//
//
#ifndef DONUTENGINE_H
#define DONUTENGINE_H

#include <string>
#include <vector>
#include <map>

#include "ArrayTypes.h"  
#include "FFTWClass.h"
#include "Zernike.h"
#include "fitsio.h"

// typedefs for DonutEngine
typedef std::vector<std::string> VString;
typedef std::vector<Complex> VComplex;
typedef std::map<std::string, std::string> MapStoS;
typedef std::map<std::string, int> MapStoI;
typedef std::map<std::string, double> MapStoD;

class DonutEngine{

public:
  // constructors
  DonutEngine(MapStoS inputMapS, MapStoI inputMapI, MapStoD inputMapD);

  // DonutEngine(int iTelescope=0,  
  // 	      double waveLength=700.0e-9,  
  // 	      int nZernikeTerms=11,  
  // 	      int nbin=256,  
  // 	      int nPixels=64,  
  // 	      const char* outputPrefix="test", 
  // 	      bool debugFlag=false, 
  // 	      int printLevel=0, 
  // 	      bool gridCalcMode=true, 
  // 	      int pixelOverSample=4, 
  // 	      double scaleFactor=2.0, 
  // 	      const char* inputPupilMask="", 
  // 	      int zemaxToDECamSignFlip=-1, 
  // 	      bool calcRzeroDerivative=false); 

  ~DonutEngine();

  // public methods - version with input arrays
  void calcAll(Real* par);
  void calcDerivatives(Real* image, Real* weight);
  void calcWFMtoImage(Matrix& wfm);
  void calcPupilFuncFromWFM(Matrix& wfm);
  void resetTimers();
  void printTimers();
  Vector& getvParCurrent(){return _parCurrent;};  
  Vector& getvDerivatives(){return _dChi2dpar;};  
  void savePar();
  void printOptions();
  void closeFits();

  // public methods - version for SWIG using numpy arrays or lists
  void calcAll(double* par, int n);
  void calcDerivatives(double* image, int ny, int nx, double* weight, int my, int mx);
  void getParCurrent(double** ARGOUTVIEW_ARRAY1, int* DIM1);  
  void getDerivatives(double** ARGOUTVIEW_ARRAY1, int* DIM1);  
  void setXYDECam(double x, double y){_xDECam = x; _yDECam = y;};
  void fillPar(double* par, int n);
  void calcWFMtoImage(double* IN_ARRAY2, int DIM1, int DIM2);

  // getter methods returning Matrix (caution!: these return a reference, so don't their object)
  Matrix& getXaxis(){return _xaxis;}
  Matrix& getYaxis(){return _yaxis;}
  Matrix& getRho(){return _rho;}
  Matrix& getTheta(){return _theta;}
  MatrixC& getPixelBox(){return _pixelBox;}
  MatrixC& getFtsPixel(){return _ftsPixel;}
  Matrix& getPupilMask(){return _pupilMask;}
  Matrix& getPupilWaveZernike(){return _pupilWaveZernike;}
  MatrixC& getPupilFunc(){return _pupilFunc;}
  Matrix& getPsfOptics(){return _psfOptics;}
  MatrixC& getFtsOptics(){return _ftsOptics;}
  Matrix& getPsfAtmos(){return _psfAtmos;}
  MatrixC& getFtsAtmos(){return _ftsAtmos;}
  Matrix& getConvOpticsAtmosPixel(){return _convOpticsAtmosPixel;}
  Matrix& getValPixelCenters(){return _valPixelCenters;}
  Matrix& getImage(){return _calcImage;}
  Zernike* getZernikeObject(){return _zernikeObject;}

  // getter methods for use with SWIG and numpy.i
  void getvXaxis(double** ARGOUTVIEW_ARRAY2, int* DIM1, int* DIM2);
  void getvYaxis(double** ARGOUTVIEW_ARRAY2, int* DIM1, int* DIM2);
  void getvRho(double** ARGOUTVIEW_ARRAY2, int* DIM1, int* DIM2);
  void getvTheta(double** ARGOUTVIEW_ARRAY2, int* DIM1, int* DIM2);
  void getvPixelBox(Complex** ARGOUTVIEW_ARRAY2, int* DIM1, int* DIM2);
  void getvFtsPixel(Complex** ARGOUTVIEW_ARRAY2, int* DIM1, int* DIM2);
  void getvPupilMask(double** ARGOUTVIEW_ARRAY2, int* DIM1, int* DIM2);
  void getvPupilWaveZernike(double** ARGOUTVIEW_ARRAY2, int* DIM1, int* DIM2);
  void getvPupilFunc(Complex** ARGOUTVIEW_ARRAY2, int* DIM1, int* DIM2);
  void getvPsfOptics(double** ARGOUTVIEW_ARRAY2, int* DIM1, int* DIM2);
  void getvFtsOptics(Complex** ARGOUTVIEW_ARRAY2, int* DIM1, int* DIM2);
  void getvPsfAtmos(double** ARGOUTVIEW_ARRAY2, int* DIM1, int* DIM2);
  void getvFtsAtmos(Complex** ARGOUTVIEW_ARRAY2, int* DIM1, int* DIM2);
  void getvConvOpticsAtmosPixel(double** ARGOUTVIEW_ARRAY2, int* DIM1, int* DIM2);
  void getvValPixelCenters(double** ARGOUTVIEW_ARRAY2, int* DIM1, int* DIM2);
  void getvImage(double** ARGOUTVIEW_ARRAY2, int* DIM1, int* DIM2);

  // don't strictly need these now
  // VString getVNames(){return parNames;}
  // VString getVTitles(){return parTitles;}
  double getScaleFactor(){return (double) _scaleFactor;};

  // set methods
  void setCalcRzeroDerivativeTrue(){_calcRzeroDerivative=true;};
  void setCalcRzeroDerivativeFalse(){_calcRzeroDerivative=false;};

  // Public variables (make some of the input variables Public
  int ipar_nEle,ipar_rzero,ipar_bkgd,ipar_ZernikeFirst,ipar_ZernikeLast;
  int npar;
  int _nbin;
  int _nPixels;
  int _pixelOverSample;
  VString parNames;
  VString parTitles;
  // calling statistics
  int nCallsCalcAll,nCallsCalcDerivative;
  int nZernikeSize;
  

protected:
  

  // methods for constructor initialization
  void init();
  void fillOptions(MapStoS inputMapS, MapStoI inputMapI, MapStoD inputMapD);

  // methods for internal use in initialization
  void calcParameters(int iT);
  void makePupilArrays(int nbin,Real lo,Real hi,Real radiusOuter);
  void makeXPsf(int nbin,Real lambdaz);
  void setupArrays();
  void setupStuff();
  void initStateMachine();
  void defineParams();

  // internal methods
  void fillPar(double* par);
  void calcPupilMask();
  void calcPupilFunc();
  void calcOptics();
  void calcAtmos();
  void calcFTPixel();
  void calcConvolute();
  void calcPixelate();

  // used for alignment of arrays
  size_t alignR;
  size_t alignC;

  // utility methods
  void itricksMGrid(int n,Real lo, Real hi, Matrix& x, Matrix& y);
  void realToComplex(Matrix& in, MatrixC& out);
  void absComplex(MatrixC& in,Matrix& out);
  void toFits(fitsfile* fptr, Matrix& in);
  void toFits(fitsfile* fptr, MatrixC& in);


  // constants
  Real _M_PI;
  Real _M_PI_4;

  // input maps
  MapStoS _inputMapS;
  MapStoI _inputMapI;
  MapStoD _inputMapD;
  
  // input parameters  
  int _iTelescope;
  Real _waveLength;
  int _nZernikeTerms;
  std::string _outputPrefix;
  bool _debugFlag;
  int _printLevel;
  bool _gridCalcMode;
  Real _scaleFactor;
  std::string _inputPupilMask;
  int _zemaxToDECamSignFlip;
  bool _calcRzeroDerivative;

  // telescope parameters, set from _iTelescope
  Real _outerRadius;
  Real _innerRadius;
  Real _zLength;
  Real _lambdaz;
  Real _fLength;
  Real _pixelSize;

  // pupil parameters (these are set by setXYDECam)
  Real _xDECam;
  Real _yDECam;

  // grid parameters
  Real _ngridperPixel;
  Real _deltaX;
  Real _pupilscale;
  Real _Lu;
  Real _Lf;
  int _nhalfPixels;

  // cpu time
  Real _timePupilMask,_timePupilFunc,_timeOptics,_timeAtmos,_timeConvolute,_timePixelate,_timeDerivatives0,_timeDerivatives1,_timeDerivatives2;


  // parameters, for internal use
  Real _nEle;
  Real _bkgd;
  Real _rzero;
  Vector _ZernikeArr;  
  Vector _parCurrent;  

  Real _last_nEle;
  Real _last_bkgd;
  Real _last_rzero;
  Vector _last_ZernikeArr; 
  Real _last_xDECam;
  Real _last_yDECam;

  // FFT arrays and plans
  MatrixC _fftInputArray,_fftOutputArray;
  MatrixC _ifftInputArray,_ifftOutputArray;
  Matrix _fftrtcInputArray;
  MatrixC _fftrtcTempArray;
  MatrixC _fftrtcOutputArray;

  fftw2dctc *_fft2PlanC;
  fftw2dctc *_ifft2PlanC;
  fftw2drtc *_fft2rtcPlanC;

  // Zernike object
  Zernike* _zernikeObject;

  // FFT grid arrays
  Matrix _xaxis,_yaxis;
  Matrix _rho,_theta; 

  // arrays for pupil
  MatrixC _pupilFunc,_pupilFuncStar;

  // arrays for Focal plane grid
  Matrix _xpsf,_ypsf;

  // Wavefront arrays and pupil function
  Matrix _pupilMask;
  Real _pupilSNorm;
  Matrix _pupilWaveZernike;
  
  // atmosphere arrays
  Matrix _rAtmos,_shftrAtmos;
  Matrix _psfAtmos;
  MatrixC _ftsAtmos;

  // Psf arrays
  MatrixC _calcG,_calcGstar;
  Matrix _psfOptics;
  MatrixC _ftsOptics;  

  // pixelization kernel arrays
  MatrixC _pixelBox;
  MatrixC _ftsPixel;

  // convolution array
  Matrix _convOpticsAtmosPixel;

  // State Machine flags
  bool _first,_statePupilMask,_statePupilFunc,_stateAtmos,_stateDerivatives;
  int _stateCounter;

  // pixel Array
  Matrix _valPixelCenters;
  Matrix _calcImage;

  // Derivatives
  Vector _dChi2dpar;

  // Fits file pointer for debug output
  fitsfile *_fptr;

  // DonutEngine for use in rzero Derivative calculation
  DonutEngine* _anotherDonutEngine;
    

};
#endif
