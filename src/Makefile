# Makefile in src

include $(P4EST_DIR)/etc/Makefile.p4est.mk

FC = gfortran
FFLAGS = -O0 -g -Wno-unused-dummy-argument -Wno-unused-variable

FORTRAN_SRC = $(FORESTCLAW)/src/fortran_source2d
FORTRAN_INC = $(FORTRAN_SRC)

fclawobjs = \
	ClawPatch.o \
	clawpack_fort.o \
	amrinit.o \
	amrrun.o \
	amrreset.o \
	subcycle_manager.o \
	amr_advance.o \
	amr_regrid.o \
	amr_output.o \
	amr_diagnostics.o \
	amr_exchangecoarse.o \
	amr_levelexchange.o \
	amr_mol.o \
	amr_single_step.o \
	amr_neighbors.o \
	amr_physbc.o \
	amr_utils.o \
	amr_options.o \
	regrid_functions.o \
	output_functions.o \
	fclaw2d_convenience.o \
	fclaw2d_solvers.o \
	fclaw2d_map.o \
        fclaw2d_p4est.o \
	forestclaw2d.o \
	$(FORTRAN_SRC)/philim.o \
	$(FORTRAN_SRC)/claw_out2.o \
	$(FORTRAN_SRC)/mapc2p.o \
	$(FORTRAN_SRC)/mapc2m.o \
	$(FORTRAN_SRC)/setup_mesh.o \
	$(FORTRAN_SRC)/claw2d_utils_mapped.o \
	$(FORTRAN_SRC)/claw2d_utils.o \
	$(FORTRAN_SRC)/debug.o \
	$(FORTRAN_SRC)/cellave2.o \
	$(FORESTCLAW)/src/mappings/set_maptype.o \
	$(FORESTCLAW)/src/mappings/scale_rot_map.o \
	$(FORESTCLAW)/src/mappings/block.o \
	$(FORESTCLAW)/src/mappings/cart/mapc2m_cart.o \
	$(FORESTCLAW)/src/mappings/disk/mapc2m_disk.o \
	$(FORESTCLAW)/src/mappings/sphere/mapc2m_sphere.o

libforestclaw.a: $(fclawobjs)
	ar cru $@ $^

.c.o:
	$(P4EST_CC) $(P4EST_CFLAGS) $(P4EST_CPPFLAGS) -c $< -o $@

.cpp.o:
	$(P4EST_CXX) $(P4EST_CFLAGS) $(P4EST_CPPFLAGS) -I. -c $< -o $@

.f.o:
	$(FC) -Wall $(FFLAGS) -c $< -o $@

clean:
	rm -f *.o *.a $(fclawobjs)

realclean:
	rm -f *.o *.a $(fclawobjs)
