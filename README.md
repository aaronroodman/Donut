# Donut
Wavefront analysis code

This code analyzes out of focus star images and retrieves the optical
wavefront, and is used for the DECam Active Optics System

It is described in Roodman etal, SPIE 2014.   Any use of this software
should cite this reference.

Setup Instructions (sorry, haven't implemented any standard build
procedures yet).

1) cd src
2) make swig
3) make or on SLAC computers: bsub -q short -o make.log make
4) add   YourArea/Donut  to PYTHONPATH
