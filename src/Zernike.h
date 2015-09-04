//
// Zernike.h:  Class to calculate Zernike polynomials of desired order
//
// Copyright (C) 2011 Aaron J. Roodman, SLAC National Accelerator Laboratory, Stanford University
//
#ifndef ZERNIKE_HH
#define ZERNIKE_HH

#include "ArrayTypes.h"

class Zernike{

public:
  // constructor
  Zernike(Matrix& rhoArr, Matrix& thetaArr,int nTerms);  
  // destructor
  virtual ~Zernike() {delete _zernikeDescription;}

  // methods
  void init(Matrix& rhoArr, Matrix& thetaArr, int nTerms);
  int factorial(int ix);

  // Public variables  
  AofMatrix _zernikeTerm;  
  std::string* _zernikeDescription;


protected:

  // input parameters  
  int _nTerms;

};
#endif
