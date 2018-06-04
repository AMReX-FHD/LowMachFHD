NDEBUG    := t
MPI       := 
OMP       :=
MKVERBOSE := t
COMP      := gfortran

# need this to compile bl_rng.f90 (module to store random number engines and distributions)
CXX11 = t

# FBOXLIB_HOME defines the directory in which we will find all the amrex code
# If you set FBOXLIB_HOME as an environment variable, this line will be ignored
FBOXLIB_HOME := ../../../../../FBoxLib
HYDROLIB_HOME := ../../../../HydroGrid
STAGGERED_HOME := ../../..

include $(FBOXLIB_HOME)/Tools/F_mk/GMakedefs.mak

include ../GPackage.mak
VPATH_LOCATIONS += ..

include $(STAGGERED_HOME)/src_gmres/GPackage.mak
VPATH_LOCATIONS += $(STAGGERED_HOME)/src_gmres

include $(STAGGERED_HOME)/src_common/GPackage.mak
VPATH_LOCATIONS += $(STAGGERED_HOME)/src_common

include $(HYDROLIB_HOME)/GPackage.mak
VPATH_LOCATIONS += $(HYDROLIB_HOME)

include $(FBOXLIB_HOME)/Src/BaseLib/GPackage.mak
VPATH_LOCATIONS += $(FBOXLIB_HOME)/Src/BaseLib

include $(FBOXLIB_HOME)/Src/MultiGrid/GPackage.mak
VPATH_LOCATIONS += $(FBOXLIB_HOME)/Src/MultiGrid

main.$(suf).exe: $(objects)
	$(LINK.f90) -o main.$(suf).exe $(objects) $(libraries)

include $(FBOXLIB_HOME)/Tools/F_mk/GMakerules.mak


