include $(HYDROLIB_HOME)/GPackageCommon.mak

myf90sources += ParallelRNG.f90

libraries += -L/usr/lib64/ -lfftw3 -lfftw3f -lm # Dynamic FFTW
#libraries += ${HPC_LIBS_64}/static/libfftw3*.a -lm # Static (manual)
#libraries += -Wl,-Bstatic -L/usr/lib64/ -lfftw3 -lfftw3f -Wl,-Bdynamic -lm # Static FFTW (using ld)
