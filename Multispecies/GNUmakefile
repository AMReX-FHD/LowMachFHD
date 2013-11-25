# BOXLIB_HOME defines the directory in which we will find all the BoxLib code
# If you set BOXLIB_HOME as an environment variable, this line will be ignored
BOXLIB_HOME = ../../../BoxLib
STAGGERED_HOME := ..

NDEBUG    := t 
MPI       :=  
OMP       :=
PROF      :=
COMP      := gfortran
MKVERBOSE := t

include $(BOXLIB_HOME)/Tools/F_mk/GMakedefs.mak

include $(STAGGERED_HOME)/Multispecies/GPackage.mak
VPATH_LOCATIONS += $(STAGGERED_HOME)/Multispecies

include $(STAGGERED_HOME)/src_common/GPackage.mak
VPATH_LOCATIONS += $(STAGGERED_HOME)/src_common

include $(BOXLIB_HOME)/Src/F_BaseLib/GPackage.mak
VPATH_LOCATIONS += $(BOXLIB_HOME)/Src/F_BaseLib

LAPACK95_DIR = $(STAGGERED_HOME)/Multispecies/LAPACK95

F90FLAGS += -I$(LAPACK95_DIR)/lapack95_modules
libraries += $(LAPACK95_DIR)/lapack95.a -llapack -lblas

main.$(suf).exe: $(objects) 
	$(LINK.f90) -o main.$(suf).exe $(objects) $(libraries)

include $(BOXLIB_HOME)/Tools/F_mk/GMakerules.mak
