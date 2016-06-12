# $Id: Makefile,v 1.1.1.1 2001-02-05 15:55:51 cbilderback Exp $
# Makefile for the UGM program
# Last modified by $Author: cbilderback $ on $Date: 2001-02-05 15:55:51 $

prefix       = /home/cbilderback
host_os      = linux
srcdir       = .
top_srcdir   = .
enable_debug = no

# Set up the include paths
INCPATHS = -I$(prefix)/Include
LIBDIRS  = -L$(prefix)/Lib

# Libraries we need to link in
LIBS = -lProjection -lgeolib -lProjectionIO -lProjectionMesh -lMathLib 

# Linker flags
LDFLAGS   = $(LIBDIRS)
LOADLIBES = $(LIBS)

# Set our compiler options
ifeq ($(enable_debug),yes)
DEBUG = -g -Wall -DTNT_NO_BOUNDS_CHECK
else
DEBUG = -O2 -DTNT_NO_BOUNDS_CHECK
#-march=pentiumpro -mcpu=pentiumpro -fomit-frame-pointer -mieee-fp -fschedule-insns2 -finline-functions -frerun-loop-opt -fstrength-reduce -ffast-math -funroll-loops -fexpensive-optimizations -fthread-jumps
endif

# Compiler and other defs
CC   = gcc
CXX  = g++ -fpermissive
CXXFLAGS = $(DEBUG) $(INCPATHS)

# Suffix rules
.SUFFIXES: .o .cpp
.cpp.o:
	$(CXX) $(CXXFLAGS) -c $<

# Dependencies for the program
OBJS = pointprojector.o

test : pointprojector.o
	$(CXX) $(CXXFLAGS) pointprojector.o -o pproject $(LIBDIRS) $(LIBS)


clean:
	rm -f $(OBJS) *~ test? test











