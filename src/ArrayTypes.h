//
// $Rev:: 191                                                       $:  
// $Author:: roodman                                                $:  
// $LastChangedDate:: 2014-09-03 11:00:33 -0700 (Wed, 03 Sep 2014)  $:  
//
// ArrayTypes.h :  defines typedefs for multidimensional arrays used in 
// Copyright (C) 2011 Aaron J. Roodman, SLAC National Accelerator Laboratory
//

#ifndef ARRAYTYPES_HH
#define ARRAYTYPES_HH

#include <complex>
#include <fftw3.h>
#include "Array.h"

typedef double Real;
typedef std::complex<double> Complex;
typedef Array::array1<Real> Vector;
typedef Array::array2<Real> Matrix;
typedef Array::array2<Complex> MatrixC;
typedef Array::array3<Real> AofMatrix;

#endif
