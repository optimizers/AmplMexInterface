include Makefile.user

CC       = gcc
MEX      = $(MATLAB)/bin/mex
MCFLAGS += -largeArrayDims -lmwblas
MEX_INTERFACES = ampl_interface_mex spam_interface_mex

AMPL    = $(shell brew --prefix asl)

all: ${MEX_INTERFACES}

${MEX_INTERFACES}:
	$(MEX) $(MCFLAGS) -I$(AMPL)/include $@.cpp -L$(AMPL)/lib -lasl -lfuncadd0

clean:
	rm -rf *.o *~

veryclean:
	rm -rf *.o *~ *.a *.mex*






