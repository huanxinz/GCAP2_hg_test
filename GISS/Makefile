#------------------------------------------------------------------------------
#                  GEOS-Chem Global Chemical Transport Model                  #
#------------------------------------------------------------------------------
#BOP
#
# !MODULE: Makefile 
#
# !DESCRIPTION: This makefile compiles the various ModelE routines necessary
#  for tracer advection for GEOS-Chem
#
# !REVISION HISTORY: 
#  10 Oct 2014 - L. Murray - Initial Version
#EOP
#------------------------------------------------------------------------------
#BOC

###############################################################################
###                                                                         ###
###  Initialization section                                                 ###
###                                                                         ###
###############################################################################

# Define variables
ROOTDIR := ..
HDR     := $(ROOTDIR)/include
LIB     := $(ROOTDIR)/lib
MOD     := $(ROOTDIR)/mod

# Include header file.  This returns CC, F90, FREEFORM, LD, R8, SHELL,
# as well as the default Makefile compilation rules for source code files.
include Makefile_header.mk

# List files needed for QUS
OBJ := giss_gc_bridge_mod.o
OBJ += QUS3D.o QUS_DRV.o QUSDEF.o QUS_COM.o MOMEN2ND.o ATMDYN.o

ifeq ($(GRID),F40)
OBJ += RES_F40.o FFT144.o 
endif

ifeq ($(GRID),M23)
OBJ += RES_M23.o FFT72.o 
endif

###############################################################################
###                                                                         ###
###  Makefile targets: type "make help" for a complete listing!             ###
###                                                                         ###
###############################################################################

.PHONY: clean help

lib: $(OBJ)
	$(AR) crs libGISS.a $(OBJ)
	mv libGISS.a $(LIB)

clean:
	rm -f *.i *.i90 *.o *.mod *.smod

help:
	@$(MAKE) -C $(HELP)

###############################################################################
###                                                                         ###
###  Dependencies listing                                                   ###
###  (grep "USE " to get the list of module references!)                    ###
###                                                                         ###
###  From this list of dependencies, the "make" utility will figure out     ###
###  correct order of compilation (so we don't have to do that ourselves).  ###
###  This also allows us to compile on multiple processors with "make -j".  ###
###                                                                         ###
###  NOTES:                                                                 ###
###  (1) Only specify object-file dependencies that are within this         ###
###       directory.  Object files in other directories will be referenced  ### 
###       at link-time.                                                     ###
###  (2) For "make -jN" (i.e. compile N files simultaneously), all files    ###
###       in this directory must have a listed dependency.                  ###
###                                                                         ###
###############################################################################

giss_gc_bridge_mod.o: resolution.mod
constant.mod: giss_gc_bridge_mod.o
dd2d_utils.mod: giss_gc_bridge_mod.o
dist_grid_mod.mod: giss_gc_bridge_mod.o
domain_decomp_1d.mod: giss_gc_bridge_mod.o
dynamics.mod: giss_gc_bridge_mod.o
gatherscatter_mod.mod: giss_gc_bridge_mod.o
geom.mod: giss_gc_bridge_mod.o
globalsum_mod.mod: giss_gc_bridge_mod.o
halo_mod.mod: giss_gc_bridge_mod.o
hidden_mod.mod: giss_gc_bridge_mod.o
model_com.mod: giss_gc_bridge_mod.o
param.mod: giss_gc_bridge_mod.o
pblcom.mod: giss_gc_bridge_mod.o
rad_cosz0.mod: giss_gc_bridge_mod.o
zonalmean_mod.mod: giss_gc_bridge_mod.o

QUS3D.o: domain_decomp_1d.mod dynamics.mod geom.mod model_com.mod quscom.mod
QUS3D.o: qusdef.mod
tracer_adv.mod: QUS3D.o

QUS_COM.o: domain_decomp_1d.mod model_com.mod
QUS_COM.o: qusdef.mod
somtq_com.mod: QUS_COM.o

QUS_DRV.o: domain_decomp_1d.mod qusdef.mod
quscom.mod: QUS_DRV.o

QUSDEF.o: domain_decomp_1d.mod
qusdef.mod: QUSDEF.o

MOMEN2ND.o: domain_decomp_1d.mod dynamics.mod geom.mod model_com.mod
moments.mod: MOMEN2ND.o

ATMDYN.o: rundeck_opts.h
ATMDYN.o: constant.mod domain_decomp_1d.mod dynamics.mod
ATMDYN.o: geom.mod model_com.mod
ATMDYN.o: moments.mod pblcom.mod resolution.mod somtq_com.mod tracer_adv.mod
#ATMDYN.o: tracer_com.mod
atmdyn.mod: ATMDYN.o
unrdrag_com.mod: ATMDYN.o

ifeq ($(GRID),F40)
FFT144.o: constant.mod
fft144.mod: FFT144.o
resolution.mod: RES_F40.o
endif

ifeq ($(GRID),M23)
FFT72.o: constant.mod
fft72.mod: FFT72.o
resolution.mod: RES_M23.o
endif

#ifeq ($(HDF5),yes)
#	$(F90) -DUSE_HDF5 -I$(HDF_INC) -c $<
#endif

#EOC
