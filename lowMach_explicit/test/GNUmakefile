NDEBUG := 
MPI    := t
OMP    :=
MKVERBOSE := t

COMP := gfortran

BOXLIB_HOME := ../../../BoxLib
HYDRO_DIR   := ../../HydroGrid

include $(BOXLIB_HOME)/Tools/F_mk/GMakedefs.mak

include ../GPackage.mak
VPATH_LOCATIONS += ..

include $(HYDRO_DIR)/GPackage.mak
VPATH_LOCATIONS += $(HYDRO_DIR)

include $(BOXLIB_HOME)/Src/F_BaseLib/GPackage.mak
VPATH_LOCATIONS += $(BOXLIB_HOME)/Src/F_BaseLib

include $(BOXLIB_HOME)/Src/LinearSolvers/F_MG/GPackage.mak
VPATH_LOCATIONS += $(BOXLIB_HOME)/Src/LinearSolvers/F_MG

main.$(suf).exe: $(objects)
	$(LINK.f90) -o main.$(suf).exe $(objects) $(libraries)

include $(BOXLIB_HOME)/Tools/F_mk/GMakerules.mak
