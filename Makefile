CC=g++

LDFLAGS=-lfftw3 -lm
CXXFLAGS=-Ofast -DNDEBUG 
CXXFLAGS+=-g -Wall -ansi -fomit-frame-pointer -fstrict-aliasing -ffast-math -msse2 -mfpmath=sse
CXXFLAGS+=-I$(FFTWPP_INCLUDE_PATH)

# Stick to one thread for now:
CXXFLAGS+=-DFFTWPP_SINGLE_THREAD


VPATH=.:$(FFTWPP_INCLUDE_PATH)

fftw++.o : fftw++.cc fftw++.h
	$(CXX) $(CXXFLAGS) -o $@ -c $< $(LDFLAGS)

timestepper.o : timestepper.cpp timestepper.hpp
	$(CXX) $(CXXFLAGS) -o $@ -c $< $(LDFLAGS)

vlafoo.o : vlafoo.cpp vlafoo.hpp
	$(CXX) $(CXXFLAGS) -o $@ -c $< $(LDFLAGS)

vlafoo: timestepper.o vlafoo.o fftw++.o
	$(CC) $(CXXFLAGS) $^ $(LDFLAGS) -o $@

make clean:
	rm -f vlafoo *.o
