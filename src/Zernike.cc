//
// Zernike.cc:  Class to calculate Zernike polynomials of desired order
//
// Copyright (C) 2011 Aaron J. Roodman, SLAC National Accelerator Laboratory, Stanford University
//
#include <iostream>
#include <cmath>
#include "Zernike.h"


// constructor 
Zernike::Zernike(Matrix& rhoArr,Matrix& thetaArr,int nTerms){
  init(rhoArr,thetaArr,nTerms);
}

void Zernike::init(Matrix& rho,Matrix& theta,int nTerms){
  // initialize Zernike class here, and fill Zernike array

  _nTerms = nTerms;

  // must have at least 3 terms! or my simple algorithm is unhappy
  if (nTerms<3) {
    nTerms = 3;
  } 

  // array dimensions - rho and theta must be the same
  int nx = rho.Nx();
  int ny = rho.Ny();

  // new member variables
  _zernikeDescription = new std::string[nTerms];
  _zernikeTerm.Dimension(nTerms,ny,nx);
  _zernikeTerm.Activate();

  // fill the text descriptions
  std::string descriptions[37] =       {"Piston (Bias)          1",
                                   "Tilt X                 4^(1/2) (p) * COS (A)",
                                   "Tilt Y                 4^(1/2) (p) * SIN (A)",
                                   "Power (Defocus)        3^(1/2) (2p^2 - 1)",
                                   "Astigmatism Y          6^(1/2) (p^2) * SIN (2A)",
                                   "Astigmatism X          6^(1/2) (p^2) * COS (2A)",
                                   "Coma Y                 8^(1/2) (3p^3 - 2p) * SIN (A)",
                                   "Coma X                 8^(1/2) (3p^3 - 2p) * COS (A)",
                                   "Trefoil Y              8^(1/2) (p^3) * SIN (3A)",
                                   "Trefoil X              8^(1/2) (p^3) * COS (3A)",
                                   "Primary Spherical      5^(1/2) (6p^4 - 6p^2 + 1)",
                                   "Secondary Astig X      10^(1/2) (4p^4 - 3p^2) * COS (2A)",
                                   "Secondary Astig Y      10^(1/2) (4p^4 - 3p^2) * SIN (2A)",
                                   "TetraFoil X            10^(1/2) (p^4) * COS (4A)",
                                   "TetraFoil Y            10^(1/2) (p^4) * SIN (4A)",
                                   "Secondary Coma X       12^(1/2) (10p^5 - 12p^3 + 3p) * COS (A)",
                                   "Secondary Coma Y       12^(1/2) (10p^5 - 12p^3 + 3p) * SIN (A)",
                                   "Secondary Trefoil X    12^(1/2) (5p^5 - 4p^3) * COS (3A)",
                                   "Secondary Trefoil Y    12^(1/2) (5p^5 - 4p^3) * SIN (3A)",
                                   "Pentafoil X            12^(1/2) (p^5) * COS (5A)",
                                   "Pentafoil Y            12^(1/2) (p^5) * SIN (5A)",
                                   "Secondary Spherical    7^(1/2) (20p^6 - 30p^4 + 12p^2 - 1)",
                                   "Tertiary Astig Y       14^(1/2) (15p^6 - 20p^4 + 6p^2) * SIN (2A)",
                                   "Tertiary Astig X       14^(1/2) (15p^6 - 20p^4 + 6p^2) * COS (2A)",
                                   "Secondary Tetrafoil Y  14^(1/2) (6p^6 - 5p^4) * SIN (4A)",
                                   "Secondary Tetrafoil X  14^(1/2) (6p^6 - 5p^4) * COS (4A)",
                                   "Sextafoil Y            14^(1/2) (p^6) * SIN (6A)",
                                   "Sextafoil X            14^(1/2) (p^6) * COS (6A)",
                                   "Tertiary Coma Y        16^(1/2) (35p^7 - 60p^5 + 30p^3 - 4p) * SIN (A)",
                                   "Tertiary Coma X        16^(1/2) (35p^7 - 60p^5 + 30p^3 - 4p) * COS (A)",
                                   "Tertiary Trefoil Y     16^(1/2) (21p^7 - 30p^5 + 10p^3) * SIN (3A)",
                                   "Tertiary Trefoil X     16^(1/2) (21p^7 - 30p^5 + 10p^3) * COS (3A)",
                                   "Secondary Pentafoil Y  16^(1/2) (7p^7 - 6p^5) * SIN (5A)",
                                   "Secondary Pentafoil X  16^(1/2) (7p^7 - 6p^5) * COS (5A)",
                                   "Septafoil Y            16^(1/2) (p^7) * SIN (7A)",
                                   "Septafoil X            16^(1/2) (p^7) * COS (7A)",
                                   "Tertiary Spherical     9^(1/2) (70p^8 - 140p^6 + 90p^4 - 20p^2 + 1)"};

  for (int iZ=0;iZ < nTerms;iZ++){
    if (iZ < 37){
      _zernikeDescription[iZ] = descriptions[iZ];
    } else {
      std::ostringstream desc;
      desc << "Zernike" << iZ << std::endl;
      _zernikeDescription[iZ] = desc.str();
    }
  }

  // now calculate the terms

  // first figure out what nMax is for nTerms
  int numZ = 0;
  int nLast = -1;

  while (numZ < nTerms){ 
    nLast = nLast + 1;
    for (int m=0;m < (nLast+1);m++){
      if (((nLast-m) % 2)==0){
	if (m==0) {
	  numZ = numZ + 1;
	} else {
	  numZ = numZ + 2;
	}
      }
    }
  }

  // set range n for  rho^n
  int nRhoPowers = nLast + 1;
  int nTrigTerms = nLast + 1;

  // first fill rho^n terms, cos(n theta) and sin(n theta)

  //
  // indexing is rhoPowersArr[0] is r^0
  //                              [1] is r^1  etc...
  //
  AofMatrix rhoNArr(nRhoPowers,nx,ny);

  // set rho^0 = 1.
  Matrix temp0 = rhoNArr[0];
  temp0 = 1.0;

  // set rho^1 = rho
  Matrix temp1 = rhoNArr[1];
  temp1 = rho;
  
  // set higher powers
  for (int iPower=2; iPower < nRhoPowers; iPower++){
    Matrix tempN = rhoNArr[iPower];
    Matrix tempNminus1 = rhoNArr[iPower-1];
    for (int i=0; i<nx*ny ; i++){
      tempN(i) = tempNminus1(i) * rho(i);  
    }
  }

  // now for cos,sin terms
  AofMatrix cosNArr(nTrigTerms,nx,ny);
  AofMatrix sinNArr(nTrigTerms,nx,ny);
       
  //
  // self.cosTheta[0] = is just 1
  // self.cosTheta[1] = cos(theta)
  // self.cosTheta[2] = cos(2 theta)  etc...
  //
  cosNArr[0] = 1.0;
  sinNArr[0] = 0.0;

  for (int iTrigTerm=0; iTrigTerm<nTrigTerms; iTrigTerm++){
    Matrix aCosTerm = cosNArr[iTrigTerm];
    Matrix aSinTerm = sinNArr[iTrigTerm];
    for  (int i=0; i<nx*ny ; i++){
      aCosTerm(i) = cos( iTrigTerm * theta(i) );
      aSinTerm(i) = sin( iTrigTerm * theta(i) );
    }
  }        

  //
  // now calculate the Zernike polynomials
  //
  int iZ = -1;   // Zernike term counter
  for (int n=0; n<nRhoPowers; n++){
    for (int m=0; m < n+1 ; m++){

      if (((n-m) % 2)==0){

	int nn = (n-m)/2;
	int nm = (n+m)/2;

	// radial term
	Matrix radialTerm(nx,ny);
	radialTerm = 0.0;
	for (int s=0; s < nn+1 ; s++){

	  Real rcoeff = pow(-1.0,s) * factorial(n-s) / (  factorial(s) * factorial(nm-s) * factorial(nn-s)   );
	  Matrix rhoNTerm = rhoNArr[n-2*s];
	  for  (int i=0; i<nx*ny ; i++){
	    radialTerm(i) = radialTerm(i) + rcoeff * rhoNTerm(i);
	  }
	}

	Real coeff = sqrt(2.0*n+2.0);
	if (m==0){
	  coeff = coeff/sqrt(2.0);
	} 
	  
	// even and odd terms
	if (m==0){
	  iZ = iZ + 1;
	  if (iZ<nTerms){
	    Matrix aZernikeTerm = _zernikeTerm[iZ];
	    for  (int i=0; i<nx*ny ; i++){
	      aZernikeTerm(i) = coeff * radialTerm(i);
	    } 
	    // std::cout << "zernike #"<<iZ<< "  n,m,coeff =" << n << " " << m << " " << coeff << std::endl; 
	  }
	  
	} else{   
	  // convention in Noll is that iZ odd is sin, and iZ even is cos
	  // but Noll also has iZ starting from 1, but we start at 0, so it is reversed!
	  
	  iZ = iZ + 1;
	  if ((iZ % 2)!=0){
	    
	    if (iZ<nTerms){
	      Matrix aZernikeTerm = _zernikeTerm[iZ];
	      Matrix aTrigTerm = cosNArr[m];
	      for  (int i=0; i<nx*ny ; i++){
		aZernikeTerm(i) = coeff * radialTerm(i) * aTrigTerm(i);
	      } 
	      // std::cout << "zernike #"<<iZ<< "  n,m,coeff =" << n << " " << m << " " << coeff << std::endl; 
	      
	    }
	    
	    iZ = iZ + 1;
	    if (iZ<nTerms){
	      Matrix aZernikeTerm = _zernikeTerm[iZ];
	      Matrix aTrigTerm = sinNArr[m];
	      for  (int i=0; i<nx*ny ; i++){
		aZernikeTerm(i) = coeff * radialTerm(i) * aTrigTerm(i);
	      }
	      // std::cout << "zernike #"<<iZ<< "  n,m,coeff =" << n << " " << m  << " " << coeff << std::endl; 
	      
	    }
	    
	  } else {

	    if (iZ<nTerms){
	      Matrix aZernikeTerm = _zernikeTerm[iZ];
	      Matrix aTrigTerm = sinNArr[m];
	      for  (int i=0; i<nx*ny ; i++){
		aZernikeTerm(i) = coeff * radialTerm(i) * aTrigTerm(i);
	      }
	      // std::cout << "zernike #"<<iZ<< "  n,m,coeff =" << n << " " << m  << " " << coeff << std::endl; 
 
	    }
	      
	    iZ = iZ + 1;
	    if (iZ<nTerms){
	      Matrix aZernikeTerm = _zernikeTerm[iZ];
	      Matrix aTrigTerm = cosNArr[m];
	      for  (int i=0; i<nx*ny ; i++){
		aZernikeTerm(i) = coeff * radialTerm(i) * aTrigTerm(i);
	      }
	      // std::cout << "zernike #"<<iZ<< "  n,m,coeff =" << n << " " << m << " " << coeff << std::endl; 
 
	    }
	      
	  } //IZ even or odd

	} // m not equal 0

      } // n-m even

    } // loop over m

  } // loop over n

}

    
int Zernike::factorial(int ix){
  int answer(0);
  if (ix==1 || ix==0){
    answer = 1;
  } else {
    answer = ix * factorial(ix-1);
  }
  return answer;
}
