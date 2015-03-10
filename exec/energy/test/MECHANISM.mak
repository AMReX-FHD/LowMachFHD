# default target (make just takes the one that appears first)
ALL: main.$(suf).exe


#-----------------------------------------------------------------------------
# core BoxLib directories
BOXLIB_CORE := $(BOXLIB_HOME)/Src/F_BaseLib


#-----------------------------------------------------------------------------

ifeq ($(strip $(BOXLIB_HOME)),)
  $(error "Must define BOXLIB_HOME")
endif

#-----------------------------------------------------------------------------
# Fmpack is the list of all the GPackage.mak files that we need to
# include into the build system to define the list of source files.
#
# Fmlocs is the list of all the directories that we want to search
# for the source files in -- this is usually going to be the
# same as the list of directories containing GPackage.mak defined
# above.
#
# Fincs is the list of directories that have include files that
# we need to tell the compiler about.

# Chemistry & Transport
Fmpack += $(CHEMISTRY_DIR)/src_common/GPackage.mak
Fmlocs += $(CHEMISTRY_DIR)/src_common
Fmpack += $(CHEMISTRY_DIR)/F_Src/GPackage.mak
Fmlocs += $(CHEMISTRY_DIR)/F_Src

INCLUDE_LOCATIONS += $(CHEMISTRY_DIR)/src_common
ifeq ($(CHEMISTRY_MODEL),INERT)
  csources += inert30.c
  vpath %.c   $(VPATH_LOCATIONS) $(CHEMISTRY_DIR)/data/inert
endif

# include the necessary GPackage.mak files that define this setup
include $(Fmpack)

# vpath defines the directories to search for the source files

#  Note: GMakerules.mak will include '.' at the start of the
#  VPATH_LOCATIONS to first search in the problem directory
VPATH_LOCATIONS += $(Fmlocs)

# list of directories to put in the Fortran include path
FINCLUDE_LOCATIONS += $(Fmincs)


#-----------------------------------------------------------------------------
# define the build instructions for the executable
main.$(suf).exe: $(objects)
ifdef MKVERBOSE
	$(HPCLINK) $(LINK.f90) -o main.$(suf).exe $(objects) $(libraries)
	@echo SUCCESS
else
	@echo "Linking $@ ..."
	@$(HPCLINK) $(LINK.f90) -o main.$(suf).exe $(objects) $(libraries)
	@echo SUCCESS
endif


#-----------------------------------------------------------------------------
# runtime parameter stuff (probin.f90)

# template used by write_probin.py to build probin.f90
PROBIN_TEMPLATE := $(SMC_TOP_DIR)/probin.template

# list of the directories to search for _parameters files
PROBIN_PARAMETER_DIRS = ./

PROBIN_PARAMETER_DIRS += $(SMC_TOP_DIR)

# list of all valid _parameters files for probin
PROBIN_PARAMETERS := $(shell $(BOXLIB_HOME)/Tools/F_scripts/findparams.py $(PROBIN_PARAMETER_DIRS))

probin.f90: $(PROBIN_PARAMETERS) $(PROBIN_TEMPLATE)
	@echo " "
	@echo "${bold}WRITING probin.f90${normal}"
ifdef MKVERBOSE
	$(BOXLIB_HOME)/Tools/F_scripts/write_probin.py \
           -t $(PROBIN_TEMPLATE) -o probin.f90 -n probin \
           --pa "$(PROBIN_PARAMETERS)"
else
	@$(BOXLIB_HOME)/Tools/F_scripts/write_probin.py \
           -t $(PROBIN_TEMPLATE) -o probin.f90 -n probin \
           --pa "$(PROBIN_PARAMETERS)"
endif
	@echo " "


#-----------------------------------------------------------------------------
# build_info stuff
deppairs: build_info.f90

build_info.f90:
	@echo " "
	@echo "${bold}WRITING build_info.f90${normal}"
ifdef MKVERBOSE
	$(BOXLIB_HOME)/Tools/F_scripts/make_build_info2 \
            "$(Fmdirs)" "$(COMP)" "$(FCOMP_VERSION)" \
            "$(COMPILE.f90)" "$(COMPILE.f)" \
            "$(COMPILE.c)" "$(LINK.f90)" \
            "AUX=$(CHEMISTRY_MODEL)" \
            "GIT=$(BOXLIB_HOME)" "GIT=$(SMC_TOP_DIR)"
else
	@$(BOXLIB_HOME)/Tools/F_scripts/make_build_info2 \
            "$(Fmdirs)" "$(COMP)" "$(FCOMP_VERSION)" \
            "$(COMPILE.f90)" "$(COMPILE.f)" \
            "$(COMPILE.c)" "$(LINK.f90)" \
            "AUX=$(CHEMISTRY_MODEL)" \
            "GIT=$(BOXLIB_HOME)" "GIT=$(SMC_TOP_DIR)"
endif
	@echo " "

$(odir)/build_info.o: build_info.f90
ifdef MKVERBOSE
	$(COMPILE.f90) $(OUTPUT_OPTION) build_info.f90
	rm -f build_info.f90
else
	@echo "Building build_info.f90 ..."
	@$(COMPILE.f90) $(OUTPUT_OPTION) build_info.f90
	@rm -f build_info.f90
endif


#-----------------------------------------------------------------------------
# include the BoxLib Fortran Makefile rules
include $(BOXLIB_HOME)/Tools/F_mk/GMakerules.mak


#-----------------------------------------------------------------------------
# for debugging.  To see the value of a Makefile variable,
# e.g. Fmlocs, simply do "make print-Fmlocs".  This will
# print out the value.
print-%: ; @echo $* is $($*)


#-----------------------------------------------------------------------------
# cleaning.  Add more actions to 'clean' and 'realclean' to remove
# probin.f90 and build_info.f90 -- this is where the '::' in make comes
# in handy
clean::
	$(RM) probin.f90
	$(RM) build_info.f90
