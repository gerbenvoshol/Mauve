
if DEBUG
D_CXXFLAGS = -Wall -g -D__GNDEBUG__ 
else
endif
OPTIMIZATION = -O3 -funroll-loops -fomit-frame-pointer -ftree-vectorize 
AM_CXXFLAGS = -DCOMMAND_LINE $(OPTIMIZATION) $(D_CXXFLAGS) 
AM_LDFLAGS = @STATIC_FLAG@ 

bin_PROGRAMS = sgEvolver scoreAlignment2 mfaToMultiFiles maf2xmfa
sgEvolver_SOURCES = \
sgEvolver.cpp Mutator.cpp Alignment.cpp \
gnGapSource.cpp gnStringSource.cpp Alignment.h \
gnStringSource.h  GisSplayTree.h \
GappedIntervalSequence.h  gnGapSource.h  Mutator.h \
twister.c twister.h

LIBRARY_CL = $(DEPS_LIBS) 

sgEvolver_LDADD = $(LIBRARY_CL)
EXTRA_sgEvolver_SOURCES = getopt.c getopt.h getopt1.c
sgEvolver_DEPENDENCIES = @GETOPT_LONG_SYSTEM@

scoreAlignment2_SOURCES = scoreAlignment2.cpp
scoreAlignment2_LDADD = $(LIBRARY_CL)

mfaToMultiFiles_SOURCES = mfaToMultiFiles.cpp
mfaToMultiFiles_LDADD = $(LIBRARY_CL)

maf2xmfa_SOURCES = maf2xmfa.cpp
maf2xmfa_LDADD = $(LIBRARY_CL)

INCLUDES = @DEPS_CFLAGS@ 

dist_bin_SCRIPTS = simujobrun.pl simujobgen.pl scoregen.pl rgradientlegend.R rgradientplot.R
dist_data_DATA = simujobparams.pm

