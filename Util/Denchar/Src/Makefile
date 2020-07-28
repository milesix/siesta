# ---
# Copyright (C) 1996-2021       The SIESTA group
#  This file is distributed under the terms of the
#  GNU General Public License: see COPYING in the top directory
#  or http://www.gnu.org/copyleft/gpl.txt .
# See Docs/Contributors.txt for a list of contributors.
# ---
#
# Makefile for stand-alone Denchar
#
# The idea is to use the code in the top Src directory as much as possible.
# This is achieved by the VPATH directive below.
# Other points to note, until we switch to a better building system:
#
#  The arch.make file is supposed to be in $(OBJDIR). This is normally
#  the top Obj, but if you are using architecture-dependent build directories
#  you might want to change this. (If you do not understand this, you do not
#  need to change anything. Power users can do "make OBJDIR=Whatever".)
#
#  If your main Siesta build used an mpi compiler, you might need to
#  define an FC_SERIAL symbol in your top arch.make, to avoid linking
#  in the mpi libraries even if we explicitly undefine MPI below.
#  
#  The dependency list at the end is overly large, but harmless

.SUFFIXES:
.SUFFIXES: .f .F .o .a  .f90 .F90

TOPDIR=.
MAIN_OBJDIR=.

VPATH=$(TOPDIR)/Util/Denchar/Src:$(TOPDIR)/Src

default: denchar

override WITH_MPI=
ARCH_MAKE=$(MAIN_OBJDIR)/arch.make
include $(ARCH_MAKE)
include $(MAIN_OBJDIR)/check_for_build_mk.mk

FC_DEFAULT:=$(FC)
FC_SERIAL?=$(FC_DEFAULT)
FC := $(FC_SERIAL)         # Make it non-recursive
#

INCFLAGS:=$(NCPS_INCFLAGS) $(FDF_INCFLAGS) $(PSML_INCFLAGS) $(INCFLAGS)

DENCHAR_UNDEF  = $(DEFS_PREFIX)-UMPI
DENCHAR_UNDEF += $(DEFS_PREFIX)-UCDF
DENCHAR_UNDEF += $(DEFS_PREFIX)-USIESTA__METIS
DENCHAR_UNDEF += $(DEFS_PREFIX)-UON_DOMAIN_DECOMP

DEFS := $(DEFS) $(DENCHAR_UNDEF)
FPPFLAGS := $(FPPFLAGS) $(DENCHAR_UNDEF)

# Uncomment the following line for debugging support
#FFLAGS=$(FFLAGS_DEBUG)


SYSOBJ=$(SYS).o

# Note that machine-specific files are now in top Src directory.
SIESTA_OBJS =  bessph.o chkdim.o dismin.o \
               dot.o  io.o iodm.o memory.o memory_log.o\
	       pxf.o radfft.o m_fft_gpfa.o volcel.o units.o

SIESTA_OBJS += precision.o parallel.o debugmpi.o m_io.o alloc.o listsc.o  \
         atmparams.o atom_options.o atmfuncs.o atm_types.o \
         radial.o spher_harm.o basis_io.o basis_types.o parallelsubs.o \
         schecomm.o sparse_matrices.o\
         files.o interpolation.o chemical.o xml.o \
         domain_decom.o printmatrix.o qsort.o mmio.o pspltm1.o \
         spatial.o m_getopts.o 

# Additional class objects
SIESTA_OBJS += m_uuid.o object_debug.o class_Sparsity.o \
	class_OrbitalDistribution.o \
	class_Data1D.o class_Data2D.o \
	class_SpData1D.o class_SpData2D.o \
        class_Geometry.o \
	class_Pair_Data1D.o \
	class_Fstack_Pair_Data1D.o \
	class_Pair_Geometry_SpData2D.o \
	class_Fstack_Pair_Geometry_SpData2D.o

LOCAL_OBJS= atompla.o colinear.o denchar.o length.o \
            local_die.o local_reinit.o local_sys.o matvect.o neighb.o \
            planed.o ranger.o readpla.o \
            redata_denchar.o rhoofr.o \
	    timer_local.o wavofr.o wrout.o alloc_handlers.o

ALL_OBJS_DENCHAR= $(SIESTA_OBJS) $(LOCAL_OBJS)

denchar: $(NCPS) $(ALL_OBJS_DENCHAR)
	$(FC) -o denchar \
	       $(LDFLAGS) $(ALL_OBJS_DENCHAR) $(FDF_LIBS) $(NCPS) $(PSML_LIBS) \
               $(XMLF90_LIBS) $(METIS_LIB)

clean:
	@echo "==> Cleaning object, library, and executable files"
	rm -f denchar *.o  *.a
	rm -f *.mod
	rm -f _tmp_deps deps.list _deps.list* protomake*

PROGS:= denchar
install: $(PROGS)
	cp -p $(PROGS) $(SIESTA_INSTALL_DIRECTORY)/bin

# Dependencies
.PHONY: dep
dep:
	@sfmakedepend --depend=obj --modext=o \
		$(addprefix $(VPATH)/,$(SIESTA_OBJS:.o=.f) $(SIESTA_OBJS:.o=.f90)) \
		$(addprefix $(VPATH)/,$(SIESTA_OBJS:.o=.F) $(SIESTA_OBJS:.o=.F90)) \
		$(VPATH)/class*.T90 \
		$(LOCAL_OBJS:.o=.f90) $(LOCAL_OBJS:.o=.F90) \
		$(LOCAL_OBJS:.o=.f) $(LOCAL_OBJS:.o=.F) || true
	@sed -i -e 's/\.T90\.o:/.o:/g' Makefile

# DO NOT DELETE THIS LINE - used by make depend
atm_types.o: precision.o radial.o
atmfuncs.o: atm_types.o local_sys.o precision.o radial.o spher_harm.o
atom_options.o: local_sys.o
basis_io.o: atm_types.o atmparams.o atom_options.o basis_types.o chemical.o
basis_io.o: local_sys.o precision.o radial.o xml.o
basis_types.o: alloc.o atmparams.o local_sys.o precision.o
bessph.o: local_sys.o precision.o
chemical.o: local_sys.o parallel.o precision.o
chkdim.o: local_sys.o
class_Data1D.o: alloc.o
class_Data2D.o: alloc.o
class_Data3D.o: alloc.o
class_Fstack_Pair_Data1D.o: class_Pair_Data1D.o class_Pair_Data1D.o
class_Fstack_Pair_Geometry_SpData2D.o: class_Pair_Geometry_SpData2D.o
class_Geometry.o: alloc.o
class_Pair_Data1D.o: class_Data1D.o class_Data1D.o
class_Pair_Geometry_SpData2D.o: class_Geometry.o class_SpData2D.o
class_SpData1D.o: class_Data1D.o class_Data1D.o class_Data1D.o class_Data1D.o
class_SpData1D.o: class_Data1D.o
class_SpData1D.o: class_OrbitalDistribution.o class_Sparsity.o
class_SpData2D.o: class_Data2D.o class_Data2D.o class_Data2D.o
class_SpData2D.o: class_OrbitalDistribution.o class_Sparsity.o
class_SpData3D.o: class_OrbitalDistribution.o class_Sparsity.o
class_Sparsity.o: alloc.o
class_TriMat.o: alloc.o
debugmpi.o: parallel.o
domain_decom.o: alloc.o local_sys.o parallel.o precision.o printmatrix.o
domain_decom.o: schecomm.o sparse_matrices.o
io.o: m_io.o
iodm.o: alloc.o files.o local_sys.o parallel.o parallelsubs.o precision.o
listsc.o: alloc.o
m_io.o: local_sys.o
memory.o: memory_log.o
memory_log.o: m_io.o parallel.o precision.o
parallelsubs.o: domain_decom.o local_sys.o parallel.o precision.o spatial.o
printmatrix.o: alloc.o
radfft.o: alloc.o bessph.o m_fft_gpfa.o precision.o
radial.o: alloc.o interpolation.o precision.o xml.o
schecomm.o: alloc.o
sparse_matrices.o: alloc.o class_Fstack_Pair_Geometry_SpData2D.o
sparse_matrices.o: class_OrbitalDistribution.o class_SpData1D.o
sparse_matrices.o: class_SpData2D.o class_SpData2D.o class_Sparsity.o
sparse_matrices.o: precision.o
spatial.o: precision.o
spher_harm.o: alloc.o local_sys.o precision.o units.o
spinorbit.o: atm_types.o atmfuncs.o atmparams.o basis_types.o m_mpi_utils.o
spinorbit.o: parallel.o parallelsubs.o precision.o pseudopotential.o
state_analysis.o: atomlist.o born_charge.o flook_siesta.o m_energies.o
state_analysis.o: m_fixed.o m_forces.o m_ntm.o m_spin.o m_stress.o
state_analysis.o: m_ts_global_vars.o m_wallclock.o parallel.o
state_analysis.o: remove_intramol_pressure.o siesta_cml.o siesta_geom.o
state_analysis.o: siesta_options.o sparse_matrices.o units.o write_subs.o
state_analysis.o: zmatrix.o
state_init.o: alloc.o atomlist.o class_Data2D.o
state_init.o: class_Fstack_Pair_Geometry_SpData2D.o
state_init.o: class_Pair_Geometry_SpData2D.o class_SpData1D.o class_SpData2D.o
state_init.o: class_SpData2D.o class_Sparsity.o create_Sparsity_SC.o
state_init.o: domain_decom.o files.o hsparse.o iodm_netcdf.o iodmhs_netcdf.o
state_init.o: iotdxv.o ioxv.o kpoint_scf.o kpoint_t.o ldau_specs.o local_sys.o
state_init.o: m_chess.o m_energies.o m_eo.o m_handle_sparse.o m_mixing.o
state_init.o: m_mixing_scf.o m_mpi_utils.o m_ncdf_siesta.o m_new_dm.o m_os.o
state_init.o: m_pivot_methods.o m_rmaxh.o m_sparse.o m_sparsity_handling.o
state_init.o: m_spin.o m_steps.o m_supercell.o m_test_io.o m_ts_electype.o
state_init.o: m_ts_global_vars.o m_ts_io.o m_ts_options.o m_ts_sparse.o
state_init.o: m_ts_tri_init.o normalize_dm.o overlap.o overlap_gradient.o
state_init.o: parallel.o proximity_check.o restructSpData2D.o siesta_cml.o
state_init.o: siesta_dicts.o siesta_geom.o siesta_options.o sparse_matrices.o
state_init.o: ts_dq.o ts_kpoint_scf.o units.o write_subs.o zmatrix.o
struct_init.o: alloc.o atmfuncs.o atomlist.o files.o iotdxv.o ioxv.o
struct_init.o: m_exp_coord.o m_iostruct.o m_mpi_utils.o m_steps.o parallel.o
struct_init.o: periodic_table.o siesta_cml.o siesta_geom.o siesta_master.o
struct_init.o: siesta_options.o units.o zmatrix.o
sys.o: parallel.o siesta_cml.o
timer.o: extrae_eventllist.o extrae_module.o local_sys.o m_timer.o parallel.o
timer.o: timer_tree.o
timer_tree.o: m_walltime.o
transition_rate.o: alloc.o fermid.o local_sys.o parallel.o parallelsubs.o
transition_rate.o: precision.o
ts_charge.o: class_OrbitalDistribution.o class_Sparsity.o geom_helper.o
ts_charge.o: m_ts_electype.o m_ts_method.o parallel.o precision.o
ts_dq.o: alloc.o class_OrbitalDistribution.o class_Sparsity.o geom_helper.o
ts_dq.o: m_interpolate.o m_ts_cctype.o m_ts_chem_pot.o m_ts_contour_eq.o
ts_dq.o: m_ts_electype.o m_ts_method.o parallel.o precision.o sorting.o
ts_dq.o: ts_charge.o units.o
ts_energies.o: class_OrbitalDistribution.o class_SpData1D.o class_SpData2D.o
ts_energies.o: class_Sparsity.o geom_helper.o m_energies.o m_spin.o
ts_energies.o: m_ts_electype.o m_ts_method.o m_ts_options.o parallel.o
ts_energies.o: precision.o sparse_matrices.o
ts_init.o: kpoint_scf.o m_fixed.o m_os.o m_ts_cctype.o m_ts_electrode.o
ts_init.o: m_ts_electype.o m_ts_gf.o m_ts_global_vars.o m_ts_method.o
ts_init.o: m_ts_options.o parallel.o siesta_options.o ts_kpoint_scf.o
ts_kpoint_scf.o: kpoint_scf.o kpoint_t.o m_spin.o m_ts_global_vars.o
ts_kpoint_scf.o: m_ts_tdir.o parallel.o precision.o siesta_options.o
ts_show_regions.o: m_region.o m_ts_electype.o m_ts_method.o parallel.o
ts_show_regions.o: precision.o units.o
typecell.o: precision.o
uncell.o: precision.o
units.o: precision.o
vacuum_level.o: moreParallelSubs.o precision.o
version.o: posix_calls.o
vmat.o: alloc.o atm_types.o atmfuncs.o atomlist.o listsc.o m_spin.o mesh.o
vmat.o: meshdscf.o meshphi.o parallel.o parallelsubs.o precision.o
vmatsp.o: alloc.o atm_types.o atmfuncs.o atomlist.o listsc.o mesh.o meshdscf.o
vmatsp.o: meshphi.o precision.o
vmb.o: local_sys.o m_fixed.o parallel.o precision.o
wavefunctions.o: atomlist.o kpoint_scf.o m_matswinvers.o m_spin.o parallel.o
wavefunctions.o: parallelsubs.o precision.o sparse_matrices.o
write_inp_wannier.o: alloc.o atmfuncs.o atomlist.o m_ntm.o m_orderbands.o
write_inp_wannier.o: mneighb.o parallel.o parallelsubs.o precision.o
write_inp_wannier.o: siesta2wannier90.o siesta_geom.o
write_md_record.o: atomlist.o m_energies.o m_kinetic.o m_steps.o m_stress.o
write_md_record.o: md_out.o parallel.o siesta_geom.o siesta_options.o units.o
write_orb_indx.o: atmfuncs.o cellsubs.o files.o precision.o
write_raw_efs.o: atmfuncs.o atomlist.o precision.o siesta_geom.o
write_subs.o: atomlist.o iofa.o m_energies.o m_forces.o m_iostruct.o m_region.o
write_subs.o: m_spin.o m_steps.o m_stress.o m_ts_electype.o m_ts_global_vars.o
write_subs.o: m_ts_method.o m_ts_options.o parallel.o precision.o siesta_cml.o
write_subs.o: siesta_geom.o siesta_options.o units.o zmatrix.o
writewave.o: alloc.o atmfuncs.o atomlist.o densematrix.o diag.o diag_option.o
writewave.o: get_kpoints_scale.o kpoint_scf.o local_sys.o m_spin.o m_spin.o
writewave.o: parallel.o parallelsubs.o precision.o siesta_geom.o units.o
xml.o: precision.o
atompla.o: local_sys.o precision.o
colinear.o: precision.o
denchar.o: basis_io.o listsc.o m_getopts.o parallel.o precision.o units.o
length.o: precision.o
local_reinit.o: files.o
matvect.o: precision.o
neighb.o: precision.o
planed.o: precision.o
ranger.o: precision.o
readpla.o: precision.o
redata_denchar.o: precision.o
rhoofr.o: atmfuncs.o chemical.o listsc.o planed.o precision.o
wavofr.o: atmfuncs.o chemical.o planed.o precision.o
wrout.o: precision.o
class_ddata1d.o: class_Data1D.o
class_ddata2d.o: class_Data2D.o
class_dspdata1d.o: class_SpData1D.o
class_dspdata2d.o: class_SpData2D.o
class_fstack_pair_ddata1d.o: class_Fstack_Pair_Data1D.o
class_fstack_pair_gdata1d.o: class_Fstack_Pair_Data1D.o
class_fstack_pair_geometry_dspdata2d.o: class_Fstack_Pair_Geometry_SpData2D.o
class_fstack_pair_sdata1d.o: class_Fstack_Pair_Data1D.o
class_gdata1d.o: class_Data1D.o
class_gspdata1d.o: class_SpData1D.o
class_idata1d.o: class_Data1D.o
class_idata2d.o: class_Data2D.o
class_ispdata1d.o: class_SpData1D.o
class_ispdata2d.o: class_SpData2D.o
class_ldata1d.o: class_Data1D.o
class_lspdata1d.o: class_SpData1D.o
class_pair_ddata1d.o: class_Pair_Data1D.o
class_pair_gdata1d.o: class_Pair_Data1D.o
class_pair_geometry_dspdata2d.o: class_Pair_Geometry_SpData2D.o
class_pair_sdata1d.o: class_Pair_Data1D.o
class_sdata1d.o: class_Data1D.o
class_sspdata1d.o: class_SpData1D.o
class_zdata1d.o: class_Data1D.o
class_zdata2d.o: class_Data2D.o
class_zspdata1d.o: class_SpData1D.o
class_zspdata2d.o: class_SpData2D.o
gpfa_core_dp.o: m_fft_gpfa.o
gpfa_core_sp.o: m_fft_gpfa.o
gpfa_fft.o: m_fft_gpfa.o
listsc_module.o: listsc.o
m_bessph.o: bessph.o
m_object_debug.o: object_debug.o
m_radfft.o: radfft.o
mtprng.o: m_uuid.o
printmat.o: printmatrix.o
sys.o: local_sys.o
