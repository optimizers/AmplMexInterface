CC      = gcc
MEX     = /Applications/MATLAB_R2014a.app/bin/mex -largeArrayDims
#MCFLAGS = -v CC=icc
#MCFLAGS = -f $(HOME)/matlab/gccR12.sh
MEX_INTERFACES = ampl_interface_mex spam_interface_mex

AMPL    = $(HOME)/local/ampl/libampl

all: ${MEX_INTERFACES}

${MEX_INTERFACES}:
	$(MEX) $(MCFLAGS) -I$(AMPL)/Src/solvers $@.cpp -L$(AMPL)/Lib -lfuncadd0 -L$(AMPL)/Lib -lampl

clean:
	rm -rf *.o *~

veryclean:
	rm -rf *.o *~ *.a *.mex*






