ARCH     = $(shell uname -m)
ARCHF    = $(shell uname -a)

## Intel compiler (2019 at least). With AMD processors cannot use -xHOST
#CXX =icpc
#CPPFLAGS  =-std=c++17 -Ofast -axSSE3 -m64 -ipo -qopenmp  $(INCLUDE)

## GCC/CLANG++ with any processor
CXX = g++
CPPFLAGS  = -std=c++17 -Ofast -flto -march=native -fopenmp  $(INCLUDE) 

ROOT = /usr/
INCLUDE = -I$(ROOT)/include
LIBS =  -lstdc++ -L$(ROOT)/lib/ -lcfitsio



OFILES = witt.o eoswrap.o cop.o main.o


.SUFFIXES: .o .f90 .cc .cpp

%.o: %.cc %.h %.hpp %.cpp
	$(CXX) $(CPPFLAGS) -c  $< -o $@

%.o: %.f90
	$(FC) $(FFLAGS) -c  $< -o $@


all: $(OFILES)
	$(CXX) $(CPPFLAGS) -o depthOpt.x *.o $(INCLUDE) $(LIBS) 


clean:
	rm -f *.o *.mod
