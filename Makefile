CC      = gcc
MEX     = /Applications/MATLAB_R2014a.app/bin/mex -largeArrayDims
#MCFLAGS = -v CC=icc
#MCFLAGS = -f $(HOME)/matlab/gccR12.sh
MEX_INTERFACES = ampl_interface_mex spam_interface_mex

AMPL    = /usr/local/Cellar/asl/20140205/

all: ${MEX_INTERFACES}

${MEX_INTERFACES}:
	$(MEX) $(MCFLAGS) -I$(AMPL)/include $@.cpp -L$(AMPL)/lib -lasl -lfuncadd0 

clean:
	rm -rf *.o *~

veryclean:
	rm -rf *.o *~ *.a *.mex*






