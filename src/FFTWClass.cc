//
// $Rev:: 191                                                       $:  
// $Author:: roodman                                                $:  
// $LastChangedDate:: 2014-09-03 11:00:33 -0700 (Wed, 03 Sep 2014)  $:  
//
// FFTWClass.cc:  C++ interface to fftw
//                (note: code borrowed from fftw++ by John Bowman)
// Copyright (C) 2011 Aaron J. Roodman, SLAC National Accelerator Laboratory
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


