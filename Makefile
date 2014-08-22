CC=g++

LDFLAGS=-lfftw3 -lm
CXXFLAGS=-O3 -DNDEBUG 
CXXFLAGS+=-g -Wall -ansi -fomit-frame-pointer -fstrict-aliasing -ffast-math -msse2 -mfpmath=sse
CXXFLAGS+=-I$(FFTWPP_INCLUDE_PATH)

# Stick to one thread for now:
CXXFLAGS+=-DFFTWPP_SINGLE_THREAD


VPATH=.:$(FFTWPP_INCLUDE_PATH)

%.o : %.cc %.h
	$(CXX) $(CXXFLAGS) -o $@ -c $< $(LDFLAGS)
%.o : %.cpp %.hpp
	$(CXX) $(CXXFLAGS) -o $@ -c $< $(LDFLAGS)

go: timestepper.o vlafou.o fftw++.o
	$(CC) $(CXXFLAGS) $^ $(LDFLAGS) -o $@

make clean:
	rm -f go *.o
