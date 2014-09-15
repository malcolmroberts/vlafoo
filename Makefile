CC=g++

LDFLAGS=
#LDFLAGS+=-L$(HOME)/boost
ifeq ($(strip $(BOOST_LIB)),)
else
LDFLAGS+=-L$(BOOST_LIB)
endif

LDFLAGS+=-lfftw3 -lm
LDFLAGS+=-lboost_program_options
CXXFLAGS=
CXXFLAGS+=-Ofast
CXXFLAGS+=-DNDEBUG 
#CXXFLAGS+=-ftree-vectorize 
#CXXFLAGS+=-ftree-vectorizer-verbose=1
CXXFLAGS+=-march=native -mtune=native
#CXXFLAGS+=-mavx
CXXFLAGS+=-g -Wall -ansi -fomit-frame-pointer -fstrict-aliasing -ffast-math -mavx -mfpmath=sse

CXXFLAGS+=-I$(FFTWPP_INCLUDE_PATH)

ifneq ($(strip $(FFTW_INCLUDE_PATH)),)
CXXFLAGS+=-I$(FFTW_INCLUDE_PATH)
endif

# Allow the specification of the Booost location:
ifeq ($(strip $(BOOST_ROOT)),)
else
CXXFLAGS+=-I$(BOOST_ROOT)
endif

# Get git branch info, define and pass to preprocessor
GITHASH=\"$(shell git rev-parse --short HEAD)\"
CXXFLAGS+=-DGITHASH=${GITHASH}
GITBRANCH=\"$(shell git rev-parse --abbrev-ref HEAD)\"
CXXFLAGS+=-DGITBRANCH=${GITBRANCH}

# Stick to one thread for now:
CXXFLAGS+=-DFFTWPP_SINGLE_THREAD

VPATH=.:$(FFTWPP_INCLUDE_PATH)

all: vlafoo 

vlafoo: timestepper.o vlafoo.o fftw++.o
	$(CC) $(CXXFLAGS) $^ $(LDFLAGS) -o $@

# The combination of .cc/.h and .cpp/.hpp means object compilation is
# hand-specified.
fftw++.o : fftw++.cc fftw++.h
	$(CXX) $(CXXFLAGS) -o $@ -c $< $(LDFLAGS)

timestepper.o : timestepper.cpp timestepper.hpp
	$(CXX) $(CXXFLAGS) -o $@ -c $< $(LDFLAGS)

vlafoo.o : vlafoo.cpp vlafoo.hpp
	$(CXX) $(CXXFLAGS) -o $@ -c $< $(LDFLAGS)

make clean:
	rm -f vlafoo *.o

FORCE:
