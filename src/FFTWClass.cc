//
// FFTWClass.cc:  C++ interface to fftw
//                (note: code extensively borrowed from fftw++ by John Bowman)
//
// Copyright (C) 1997-2010 John C. Bowman

// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation; either version 2 of the License, or
// (at your option) any later version.

// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.

// You should have received a copy of the GNU General Public License
// along with this program; if not, write to the Free Software
// Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA. */

//
// Aaron J. Roodman, SLAC National Accelerator Laboratory, Stanford University, 2011.
//

#include "FFTWClass.h"

void fftShift(Matrix& in,Matrix& out){

  // only works for even n !!!!
  // note that for even n fftshift and ifftshift are identical
  int nx = in.Nx();
  int ny = in.Ny();

  int ny2, nx2;
  Real tmp1, tmp3;

  ny2 = ny / 2;    // half of row dimension
  nx2 = nx / 2;    // half of column dimension

// interchange entries in 4 quadrants
//   start    3  4    end    2  1      
//            1  2           4  3
//
//
  int index1(0);  // remember: index = j*nx + i 
  int index2(0);
  int index3(0);
  int index4(0);

  for (int j = 0; j < ny2; j++) {
    index1 = j*nx;
    index4 = (j+ny2)*nx + nx2;

    index2 = j*nx + nx2;
    index3 = (j+ny2)*nx;

    for (int i = 0; i < nx2; i++) {
      // use a temporary so this method can be called with in=out for shift in place
      tmp1 = in(index1);
      out(index1) = in(index4);
      out(index4) = tmp1;
    
      tmp3         = in(index3);
      out(index3)  = in(index2);
      out(index2)    = tmp3;

      // go to next entry
      index1++;
      index4++;

      index2++;
      index3++;

    }
  }

}

void fftShift(Matrix& inout){
  // shift in place
  fftShift(inout,inout);
}



void fftShift(MatrixC& in,MatrixC& out){

   // only works for even n !!!!
   // note that for even n fftshift and ifftshift are identical
   int nx = in.Nx();
   int ny = in.Ny();
   
   int ny2, nx2;
   Complex tmp1, tmp3;

   ny2 = ny / 2;    // half of row dimension
   nx2 = nx / 2;    // half of column dimension
   
   // interchange entries in 4 quadrants
   //   start    3  4    end    2  1      
   //            1  2           4  3
   //
   //
   int index1(0);  // remember: index = j*nx + i 
   int index2(0);
   int index3(0);
   int index4(0);
   
   for (int j = 0; j < ny2; j++) {
     index1 = j*nx;
     index4 = (j+ny2)*nx + nx2;
     
     index2 = j*nx + nx2;
     index3 = (j+ny2)*nx;
     
     for (int i = 0; i < nx2; i++) {
       
       // use a temporary so this method can be called with in=out for shift in place
       tmp1 = in(index1);
       out(index1) = in(index4);
       out(index4) = tmp1;
       
       tmp3         = in(index3);
       out(index3)  = in(index2);
       out(index2)    = tmp3;
       
       // go to next entry
       index1++;
       index4++;
       
       index2++;
       index3++;
       
     }
   }
   
 }
 
void fftShift(MatrixC& inout){
   // shift in place
   fftShift(inout,inout);
}


