#
# Makefile for DonutEngine
#
OS := $(shell uname)

CC = gcc
CXX = g++
SWIGFLAGS = -python -c++ -w490 -I$cwd

# define this for faster operation in Array
NDEBUG = True

ifeq ($(OS), Darwin)
#	CFLAGS = -Wall -ansi -g
	CFLAGS =  -Wall -ansi -O3  -m64 -funroll-loops -fomit-frame-pointer -ffast-math -mfpmath=sse -msse2 -mtune=native -fPIC
	LD = gcc
	LDFLAGS = -bundle -flat_namespace -undefined suppress
	LIBS = -L/opt/local/lib  -lcfitsio -lfftw3 -lm
	INCS = -I/Users/roodman/Astrophysics/Code/donutlib -I/opt/local/include -I$(PYTHONHOMEDIR)/include/python2.7 -I$(PYTHONHOMEDIR)/lib/python2.7/site-packages/numpy/core/include
	SW = swig

else
	#CFLAGS = -Wall -ansi -g -fPIC
# standard	CFLAGS =  -ansi -O3  -m64 -funroll-loops -fomit-frame-pointer -ffast-math -mfpmath=sse -msse2  -fPIC
# new, faster, use in queue CFLAGS =  -ansi -O3  -m64 -march=native -mtune=native -fPIC
	#CFLAGS =  -ansi -O3  -m64 -funroll-loops -fomit-frame-pointer -ffast-math -mfpmath=sse -msse2  -fPIC	
	CFLAGS =  -ansi -O3  -m64 -march=native -mtune=native -fPIC
	LD = g++
	LDFLAGS = -shared -export-dynamic

	LIBS = -L$(XRAY_SOFTDIR)/cfitsio/3.29/lib -L$(XRAY_SOFTDIR)/fftw/3.3/lib  -L/usr/lib64 -L$(XRAY_SOFTDIR)/python/2.7.3/lib -lcfitsio -lfftw3 -lm
	INCS = -I$(XRAY_SOFTDIR)/cfitsio/3.29/include -I$(XRAY_SOFTDIR)/fftw/3.3/include -I$(XRAY_SOFTDIR)/python/2.7.3/include/python2.7  -I$(XRAY_SOFTDIR)/python/2.7.3/lib/python2.7/site-packages/numpy/core/include

	SW = $(XRAY_SOFTDIR)/swig/2.0.4/bin/swig
endif

all: donutengine

donutengine: DonutEngine.cc
	$(CXX) -c $(CFLAGS) DonutEngine.cc $(INCS) -o DonutEngine.o
	$(CXX) -c $(CFLAGS) Zernike.cc $(INCS) -o Zernike.o
	$(CXX) -c $(CFLAGS) FFTWClass.cc $(INCS) -o FFTWClass.o
	$(CXX)  $(CFLAGS) $(INCS) -c -o DonutEngineWrap.o DonutEngineWrap.cxx 
	$(LD) $(LDFLAGS) -o _donutengine.so  DonutEngineWrap.o DonutEngine.o Zernike.o FFTWClass.o  $(LIBS)

swig:
	$(SW) $(SWIGFLAGS) -o DonutEngineWrap.cxx DonutEngine.i

clean:
	rm -f *.o
	rm -f *.cxx
	rm -f *.pyc
	rm -f *.so
