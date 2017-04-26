# AMREX_HOME defines the directory in which we will find all the amrex code
# If you set AMREX_HOME as an environment variable, this line will be ignored
AMREX_HOME := ../../../../../amrex
HYDROLIB_HOME := ../../../../HydroGrid
STAGGERED_HOME := ../../..

NDEBUG    := 
MPI       := t
OMP       :=
PROF      :=
COMP      := gfortran
CCOMP     := gcc
MKVERBOSE := t

# need this to compile bl_rng.f90 (module to store random number engines and distributions)
CXX11 = t

include $(AMREX_HOME)/Tools/F_mk/GMakedefs.mak

include $(STAGGERED_HOME)/src_reactDiff/GPackage.mak
VPATH_LOCATIONS += $(STAGGERED_HOME)/src_reactDiff

include $(STAGGERED_HOME)/src_chemistry/GPackage.mak
VPATH_LOCATIONS += $(STAGGERED_HOME)/src_chemistry

include $(STAGGERED_HOME)/src_common/GPackage.mak
VPATH_LOCATIONS += $(STAGGERED_HOME)/src_common

include $(HYDROLIB_HOME)/GPackage.mak
VPATH_LOCATIONS += $(HYDROLIB_HOME)

include $(AMREX_HOME)/Src/F_BaseLib/GPackage.mak
VPATH_LOCATIONS += $(AMREX_HOME)/Src/F_BaseLib

include $(AMREX_HOME)/Src/LinearSolvers/F_MG/GPackage.mak
VPATH_LOCATIONS += $(AMREX_HOME)/Src/LinearSolvers/F_MG

main.$(suf).exe: $(objects) 
	$(LINK.f90) -o main.$(suf).exe $(objects) $(libraries)

include $(AMREX_HOME)/Tools/F_mk/GMakerules.mak
