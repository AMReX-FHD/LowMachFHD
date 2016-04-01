# BOXLIB_HOME defines the directory in which we will find all the BoxLib code
# If you set BOXLIB_HOME as an environment variable, this line will be ignored
BOXLIB_HOME := ../../../../../BoxLib
HYDROLIB_HOME := ../../../../HydroGrid
STAGGERED_HOME := ../../..

NDEBUG    :=
MPI       := t
OMP       :=
PROF      :=
COMP      := gfortran
CCOMP     := gcc
MKVERBOSE := t

include $(BOXLIB_HOME)/Tools/F_mk/GMakedefs.mak

VPATH_LOCATIONS += ..

#######################################
# choose the EOS and model here

# ideal gas EOS, hard sphere transport
include $(STAGGERED_HOME)/src_energy/ideal_HS/GPackage.mak
VPATH_LOCATIONS += $(STAGGERED_HOME)/src_energy/ideal_HS

#include $(STAGGERED_HOME)/src_energy/SRK_HS/GPackage.mak
#VPATH_LOCATIONS += $(STAGGERED_HOME)/src_energy/SRK_HS

#######################################

include $(STAGGERED_HOME)/src_energy/GPackage.mak
VPATH_LOCATIONS += $(STAGGERED_HOME)/src_energy

include $(STAGGERED_HOME)/src_multiSpec/GPackage.mak
VPATH_LOCATIONS += $(STAGGERED_HOME)/src_multiSpec

include $(STAGGERED_HOME)/src_gmres/GPackage.mak
VPATH_LOCATIONS += $(STAGGERED_HOME)/src_gmres

include $(STAGGERED_HOME)/src_common/GPackage.mak
VPATH_LOCATIONS += $(STAGGERED_HOME)/src_common

include $(HYDROLIB_HOME)/GPackage.mak
VPATH_LOCATIONS += $(HYDROLIB_HOME)

include $(BOXLIB_HOME)/Src/F_BaseLib/GPackage.mak
VPATH_LOCATIONS += $(BOXLIB_HOME)/Src/F_BaseLib

include $(BOXLIB_HOME)/Src/LinearSolvers/F_MG/GPackage.mak
VPATH_LOCATIONS += $(BOXLIB_HOME)/Src/LinearSolvers/F_MG

LAPACK95_DIR = $(STAGGERED_HOME)/src_multiSpec/LAPACK95

F90FLAGS += -I$(LAPACK95_DIR)/lapack95_modules
ifeq ($(findstring edison, $(HOST)), edison)
  libraries += $(LAPACK95_DIR)/lapack95.a
else ifeq ($(findstring cori, $(HOST)), cori)
  libraries += $(LAPACK95_DIR)/lapack95.a
else
  libraries += $(LAPACK95_DIR)/lapack95.a -llapack -lblas
endif

main.$(suf).exe: $(objects) 
	$(LINK.f90) -o main.$(suf).exe $(objects) $(libraries)

include $(BOXLIB_HOME)/Tools/F_mk/GMakerules.mak
