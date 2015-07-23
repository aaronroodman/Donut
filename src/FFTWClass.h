//
// $Rev:: 191                                                       $:  
// $Author:: roodman                                                $:  
// $LastChangedDate:: 2014-09-03 11:00:33 -0700 (Wed, 03 Sep 2014)  $:  
//
// FFTWClass.h:  C++ interface to fftw
//                (note some code borrowed from fftw++ by John Bowman)
//
// 
/* Fast Fourier transform C++ header class for the FFTW3 Library
Copyright (C) 2004-10 John C. Bowman, University of Alberta

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License as published by
the Free Software Foundation; either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA. */ 
//
//
//
// Copyright (C) 2011 Aaron J. Roodman, SLAC National Accelerator Laboratory
//
#ifndef FFTWCLASS_HH
#define FFTWCLASS_HH

#include <cstdlib>
#include <fstream>
#include <iostream>
#include <cerrno>
#include <cmath>
#include <time.h>
#include <complex>

#include <fftw3.h>
#include "ArrayTypes.h"
#include "Array.h"




void fftShift(Matrix& in, Matrix& out);
void fftShift(Matrix& in);
void fftShift(MatrixC& in, MatrixC& out);
void fftShift(MatrixC& in);


class fftw2dctc {
protected:
  int sign;
  fftw_plan plan;
  Real cputime;
  
public:
  fftw2dctc(MatrixC& in, MatrixC& out, int sign0) {
    sign = sign0;
    plan = fftw_plan_dft_2d(in.Ny(),in.Nx(),(fftw_complex *) in(), (fftw_complex *) out(),sign,FFTW_MEASURE);  
    cputime = 0.0;
  }
  
  virtual ~fftw2dctc() {
    //    if(plan) fftw_destroy_plan(plan);
    //    std::cout << " FFT Processing time is " << cputime << std::endl;
  }
  
  void execute() {
    clock_t start = clock();
    fftw_execute(plan);
    clock_t stop = clock();
    cputime += (stop-start)/(Real)CLOCKS_PER_SEC;
  }
    

};



class fftw2drtc {
protected:
  Matrix* in;
  MatrixC* temp;
  MatrixC* out;
  fftw_plan plan;
  Real cputime;
  
public:
  fftw2drtc(Matrix& in0, MatrixC& temp0, MatrixC& out0) {
    in = &in0;
    temp = &temp0;
    out = &out0;
    plan = fftw_plan_dft_r2c_2d(in0.Ny(),in0.Nx(),(double *) in0(), (fftw_complex *) temp0(),FFTW_MEASURE);  
    cputime = 0.0;
  }
  
  virtual ~fftw2drtc() {
    //    if(plan) fftw_destroy_plan(plan);
    //std::cout << " FFT Processing time is " << cputime << std::endl;
  }
  
  void execute() {
    clock_t start = clock();

    fftw_execute(plan);

    // for r->c FFTs the output c array has the correct components for iy: 0->Ny AND ix: 0->Nx/2
    // so the "out" array needs copies of the "temp" components for iy: 0->Ny AND ix: 0->Nx/2
    // and complex conjugates of the remaining components

    // if in and out are Ny by Nx, then temp must be Ny by (Nx/2 + 1)
    int Nx = in->Nx();
    int Nxhalf = (Nx/2); //probably only works for even N
    int Nxhalfp = Nxhalf+1;

    int Ny = in->Ny();
    int Nyhalf = (Ny/2); //probably only works for even N

    
    // special cases for ix=0 and iy=0
    // for iy=0, out(0,Nx-ix) = temp*(0,ix)  and out(0,ix) = temp(0,ix)
    // for ix=0, out(Ny-iy,0) = temp*(iy,0)  and out(iy,0) = temp(iy,0)
    // fill in the missing elements in the array: out

    int index(0);
    int indexStar(0);
    int indexTemp(0);


    index = 0;
    indexTemp = 0;
    indexStar = Nx;
    for (int ix=0;ix<=Nxhalf;ix++){
      (*out)(index) = (*temp)(indexTemp);
      (*out)(indexStar) = conj((*temp)(indexTemp));
      index++;
      indexTemp++;
      indexStar--;
    }
    index = Nx;
    indexTemp = Nxhalfp;
    indexStar = (Ny-1)*Nx;
    for (int iy=1;iy<=Nyhalf;iy++){
      (*out)(index) = (*temp)(indexTemp);
      (*out)(indexStar) = conj((*temp)(indexTemp));
      index += Nx;
      indexTemp += Nxhalfp;
      indexStar -= Nx;
    }

    // now fill in values for iy: 1->Ny and ix: 1->Nxhalf
    // set out(iy,ix) = temp(iy,ix);
    // set out(Ny-iy,Nx-ix) = temp*(iy,ix);

    index = Nx;
    indexTemp = Nxhalfp;
    indexStar = (Ny-1)*Nx + (Nx);
    for (int iy=1;iy<Ny;iy++){
      index++;   //since we start the next loop at x=1
      indexTemp++;
      indexStar--;
      for (int ix=1;ix<=Nxhalf;ix++){
    	(*out)(index) = (*temp)(indexTemp);
    	(*out)(indexStar) = conj((*temp)(indexTemp));
	index++;
	indexTemp++;
	indexStar--;
      }
      index += (Nx-Nxhalfp);
      indexStar -= (Nx-Nxhalfp);
    }

    clock_t stop = clock();
    cputime += (stop-start)/(Real)CLOCKS_PER_SEC;
  }
    
};

#endif
