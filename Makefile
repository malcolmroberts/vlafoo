CXX=g++

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
# depends on .git/FETCH_HEAD?
GITHASH=\"$(shell git rev-parse --short HEAD)\"
CXXFLAGS+=-DGITHASH=${GITHASH}
# depends on .git/HEAD?
GITBRANCH=\"$(shell git rev-parse --abbrev-ref HEAD)\"
CXXFLAGS+=-DGITBRANCH=${GITBRANCH}

# Stick to one thread for now:
CXXFLAGS+=-DFFTWPP_SINGLE_THREAD

VPATH=.:$(FFTWPP_INCLUDE_PATH)

# Define the source files and objects
SRCS_CPP = clopts.cpp timestepper.cpp vlafoo.cpp
SRCS_CC = fftw++.cc
OBJS=$(SRCS_CPP:.cpp=.o)
OBJS+= $(SRCS_CC:.cc=.o)

all: vlafoo 

vlafoo: $(OBJS)
	$(CXX) $(CXXFLAGS) $^ $(LDFLAGS) -o $@

# If we are not "make clean", then creat the dependency files
ifeq ($(MAKECMDGOALS),clean)
DEPS=
else
MAKEDEPEND=$(CXXFLAGS) -O0 -M -DDEPEND
DEPS=$(SRCS_CPP:.cpp=.d)
DEPS+=":"$(SRCS_CC:.cc=.d)
.cpp.d:
	echo -n "$*.d ">$*.d
	$(CXX) $(MAKEDEPEND) $*.cpp>>$*.d
.cc.d:
	echo -n "$*.d ">$*.d
	$(CXX) $(MAKEDEPEND) $*.cc>>$*.d
-include $(DEPS)
endif

# Create the objects based on the dependency fiels
.o: %.d
	@echo $@
	$(CXX) $(CXXFLAGS) $^ -c -o $@

make clean:
	rm -f vlafoo *.o *.d

# Needed so that the .d files are generated.
.SUFFIXES: .cc .cpp .o .d

FORCE:
