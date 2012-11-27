
mol_solver := $(FORESTCLAW)/src/amr_rkc.o parabolic_step_rkc.o rkc.o

amrclawobjs_orig = \
	$(FORESTCLAW)/AMRClaw/src/philim.o \
	$(FORESTCLAW)/AMRClaw/src/limiter.o \
	$(FORESTCLAW)/AMRClaw/src2d/mapc2p.o \
	$(FORESTCLAW)/AMRClaw/src2d/b4step2.o \
	$(FORESTCLAW)/AMRClaw/src2d/setaux.o \
	$(FORESTCLAW)/AMRClaw/src2d/b4step2_mapped.o \
	$(FORESTCLAW)/AMRClaw/src2d/setaux_mapped.o \
	$(FORESTCLAW)/AMRClaw/src2d/qinit_mapped.o \
	$(FORESTCLAW)/AMRClaw/src2d/src2.o \
	$(FORESTCLAW)/AMRClaw/src2d/bc2.o \
	$(FORESTCLAW)/AMRClaw/src2d/fdisc.o \
	$(FORESTCLAW)/AMRClaw/src2d/setup_mesh.o \
	$(FORESTCLAW)/AMRClaw/src2d/claw2d_utils.o \
	$(FORESTCLAW)/AMRClaw/src2d/ClawPatch2.o \
	$(FORESTCLAW)/AMRClaw/src2d/cellave2.o \
	$(FORESTCLAW)/AMRClaw/src2d/flux2.o \
	$(FORESTCLAW)/AMRClaw/src2d/qad2.o \
	$(FORESTCLAW)/AMRClaw/src2d/step2.o \
	$(FORESTCLAW)/AMRClaw/src2d/claw_out2.o \
	$(FORESTCLAW)/AMRClaw/src2d/set_maptype.o \
	$(FORESTCLAW)/AMRClaw/src2d/scale_rot_map.o \
	$(FORESTCLAW)/AMRClaw/src2d/claw2d_utils_mapped.o \
	$(FORESTCLAW)/AMRClaw/src2d/mappings/block.o \
	$(FORESTCLAW)/AMRClaw/src2d/mappings/mapc2m_cart.o \
	$(FORESTCLAW)/AMRClaw/src2d/mappings/mapc2m_disk.o \
	$(FORESTCLAW)/AMRClaw/src2d/mappings/mapc2m_mb_sphere.o
