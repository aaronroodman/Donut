%module(docstring="donutengine calculates out of focus stars from a pupil plane Zernike expansion, for the DECam, Aaron Roodman SLAC National Accelerator Laboratory, Stanford University, 2012") donutengine

// make a docstring for Swig created code
%feature("autodoc", "3");

%{
#define SWIG_FILE_WITH_INIT
#include "DonutEngine.h"
%}

// to wrap std library objects
%include "std_vector.i"
%include "std_string.i"
%include "std_map.i"

// Instantiate templates used by example
namespace std {
   %template(IntVector) vector<int>;
   %template(DoubleVector) vector<double>;
   %template(StringVector) vector<string>;
   %template(StringToStringMap) map<string, string>;
   %template(StringToIntMap) map<string, int>;
   %template(StringToDoubleMap) map<string, double>;
}


// for numpy to c++ conversions
%include "numpy.i"
%include "stl.i"
%init %{
  import_array();
%}


// for complex types with numpy
%numpy_typemaps(Complex , NPY_CDOUBLE, int)

// numpy arguments 
%apply (double* IN_ARRAY2, int DIM1, int DIM2) {(double* image, int ny, int nx)};
%apply (double* IN_ARRAY2, int DIM1, int DIM2) {(double* weight, int my, int mx)};
%apply (double* IN_ARRAY1, int DIM1) {(double* par, int n)};
%apply (double* IN_ARRAY2, int DIM1, int DIM2) {(double* wfm, int nx, int ny)};
%apply (double* IN_ARRAY1, int DIM1, double* IN_ARRAY2, int DIM2, int DIM3) {(double* par, int n, double* dwfm, int nx, int ny)};


// Include the header file to be wrapped
%include "DonutEngine.h"


/* Rewrite the high level interface to DonutEngine, but call it donutengine */
%pythoncode %{

def donutengine(**inputDict):
  """  donutengine class for calculating out-of-focus star images from Zernike pupil basis """

  # special code for scaleFactor - be sure it is a float
  if "scaleFactor" in inputDict:
    inputDict["scaleFactor"] = float(inputDict["scaleFactor"])

  # split inputDict into 3 dictionaries - S,I,D
  paramDictS = {}
  paramDictI = {}
  paramDictD = {}
  for key in inputDict.keys():
    value = inputDict[key]
    if isinstance(value,str):
      paramDictS[key] = value
    elif isinstance(value,int):
      paramDictI[key] = value
    elif isinstance(value,float):
      paramDictD[key] = value
   
  # setup the fit engine
  myDE = DonutEngine(paramDictS,paramDictI,paramDictD)

  return myDE

%}

