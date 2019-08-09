# BOXLIB_HOME defines the directory in which we will find all the BoxLib code
# If you set BOXLIB_HOME as an environment variable, this line will be ignored
BOXLIB_HOME := ../..//../BoxLib
HYDROLIB_HOME := ../

NDEBUG    := 
MPI       := t
OMP       :=
PROF      :=
COMP      := gfortran
CCOMP     := gcc
MKVERBOSE := t

# need this to compile bl_rng.f90 (module to store random number engines and distributions)
CXX11 = t

include $(BOXLIB_HOME)/Tools/F_mk/GMakedefs.mak

include ./GPackage.mak
VPATH_LOCATIONS += .

include $(HYDROLIB_HOME)/GPackage.mak
VPATH_LOCATIONS += $(HYDROLIB_HOME)

include $(BOXLIB_HOME)/Src/F_BaseLib/GPackage.mak
VPATH_LOCATIONS += $(BOXLIB_HOME)/Src/F_BaseLib

main.$(suf).exe: $(objects) 
	$(LINK.f90) -o main.$(suf).exe $(objects) $(libraries)

include $(BOXLIB_HOME)/Tools/F_mk/GMakerules.mak
