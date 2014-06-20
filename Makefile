CC      = gcc
MEX     = /Applications/MATLAB_R2014a.app/bin/mex
MCFLAGS = -largeArrayDims -lmwblas
#MCFLAGS = -f $(HOME)/matlab/gccR12.sh
MEX_INTERFACES = ampl_interface_mex spam_interface_mex

AMPL    = $(shell brew --prefix asl)

all: ${MEX_INTERFACES}

${MEX_INTERFACES}:
	$(MEX) $(MCFLAGS) -I$(AMPL)/include $@.cpp -L$(AMPL)/lib -lasl -lfuncadd0

clean:
	rm -rf *.o *~

veryclean:
	rm -rf *.o *~ *.a *.mex*






