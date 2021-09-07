/*
Copyright (c) 2012-2021 Carsten Burstedde, Donna Calhoun, Scott Aiton
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

 * Redistributions of source code must retain the above copyright notice, this
list of conditions and the following disclaimer.
 * Redistributions in binary form must reproduce the above copyright notice,
this list of conditions and the following disclaimer in the documentation
and/or other materials provided with the distribution.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED.  IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

#include <fclaw3dx_clawpatch.h>
#include <fclaw3dx_clawpatch.hpp>
#include <fclaw3dx_clawpatch_options.h>
#include <fclaw3dx_clawpatch46_fort.h>
#include <fclaw3dx_clawpatch_output_ascii.h>
#include <fclaw2d_global.h>
#include <fclaw2d_domain.h>
#include <fclaw2d_patch.h>
#include <fclaw2d_convenience.h>
#include <fclaw2d_options.h>
#include <test/catch.hpp>
#include <test/test.hpp>
#include <fstream>

#define CHECK_BOX_DIMENSIONS(abox, mbc, mx, my, mz, mfields)\
{\
	CHECK(abox.fields() == mfields);\
	CHECK(abox.box().boxDim() == 3);\
	CHECK(abox.box().smallEnd(0) == 1-mbc);\
	CHECK(abox.box().smallEnd(1) == 1-mbc);\
	CHECK(abox.box().smallEnd(2) == 1-mbc);\
	CHECK(abox.box().bigEnd(0) == mx+mbc);\
	CHECK(abox.box().bigEnd(1) == my+mbc);\
	CHECK(abox.box().bigEnd(2) == mz+mbc);\
}
#define CHECK_BOX_EMPTY(abox)\
{\
	CHECK(abox.fields() == 0);\
	CHECK(abox.box().boxDim() == 0);\
}
TEST_CASE("fclaw3dx_clawpatch_vtable_initialize","[fclaw3dx][clawpatch]")
{
	fclaw3dx_clawpatch_vtable_initialize(4);

	fclaw3dx_clawpatch_vtable_t * clawpatch_vt = fclaw3dx_clawpatch_vt();

	CHECK(clawpatch_vt->set_user_data             == NULL);

	//ghost filling
	CHECK(clawpatch_vt->fort_copy_face              == &FCLAW3DX_CLAWPATCH46_FORT_COPY_FACE);
	CHECK(clawpatch_vt->fort_average_face           == &FCLAW3DX_CLAWPATCH46_FORT_AVERAGE_FACE);
	CHECK(clawpatch_vt->fort_interpolate_face       == &FCLAW3DX_CLAWPATCH46_FORT_INTERPOLATE_FACE);
	CHECK(clawpatch_vt->fort_copy_corner            == &FCLAW3DX_CLAWPATCH46_FORT_COPY_CORNER);
	CHECK(clawpatch_vt->fort_average_corner         == &FCLAW3DX_CLAWPATCH46_FORT_AVERAGE_CORNER);
	CHECK(clawpatch_vt->fort_interpolate_corner     == &FCLAW3DX_CLAWPATCH46_FORT_INTERPOLATE_CORNER);

	//regridding
	CHECK(clawpatch_vt->fort_tag4refinement         == &FCLAW3DX_CLAWPATCH46_FORT_TAG4REFINEMENT);
	CHECK(clawpatch_vt->fort_tag4coarsening         == &FCLAW3DX_CLAWPATCH46_FORT_TAG4COARSENING);
	CHECK(clawpatch_vt->fort_user_exceeds_threshold == NULL);
	CHECK(clawpatch_vt->fort_interpolate2fine       == &FCLAW3DX_CLAWPATCH46_FORT_INTERPOLATE2FINE);
	CHECK(clawpatch_vt->fort_average2coarse         == &FCLAW3DX_CLAWPATCH46_FORT_AVERAGE2COARSE);

	//ascii output
	CHECK(clawpatch_vt->time_header_ascii           == &fclaw3dx_clawpatch_time_header_ascii);
	CHECK(clawpatch_vt->fort_header_ascii           == &FCLAW3DX_CLAWPATCH46_FORT_HEADER_ASCII);
	CHECK(clawpatch_vt->cb_output_ascii             == &fclaw3dx_clawpatch_output_ascii_cb);
	CHECK(clawpatch_vt->fort_output_ascii           == &FCLAW3DX_CLAWPATCH46_FORT_OUTPUT_ASCII);

	//time interpolation
	CHECK(clawpatch_vt->fort_timeinterp             == &FCLAW3DX_CLAWPATCH46_FORT_TIMEINTERP);

	//ghot packing
	CHECK(clawpatch_vt->fort_local_ghost_pack       == &FCLAW3DX_CLAWPATCH46_FORT_LOCAL_GHOST_PACK);
	CHECK(clawpatch_vt->local_ghost_pack_aux        == NULL);

	//diagnostics
	CHECK(clawpatch_vt->conservation_check          == NULL);
	CHECK(clawpatch_vt->compute_error               == NULL);
	CHECK(clawpatch_vt->fort_compute_patch_error    == NULL);
	CHECK(clawpatch_vt->fort_conservation_check     == NULL);
	CHECK(clawpatch_vt->fort_compute_error_norm     == NULL);
	CHECK(clawpatch_vt->fort_compute_patch_area     == NULL);

	CHECK(clawpatch_vt->is_set                      == 1);

	fclaw2d_patch_vtable_t * patch_vt = fclaw2d_patch_vt();
	//create delete build
	//TODO document patch_vt and expose these as part to public api
	CHECK(patch_vt->patch_new                      != NULL);
	CHECK(patch_vt->patch_delete                   != NULL);
	CHECK(patch_vt->build                          != NULL);
	CHECK(patch_vt->build_from_fine                != NULL);
	CHECK(patch_vt->setup                          == NULL);
}
TEST_CASE("fclaw3dx_clawpatch patch_build","[fclaw3dx][clawpatch]")
{
	fclaw2d_global_t* glob = fclaw2d_global_new();
	fclaw_options_t fopts;
	memset(&fopts, 0, sizeof(fopts));
	fopts.mi=1;
	fopts.mj=1;
	fopts.minlevel=0;
	fopts.manifold=false;
	fopts.bx = 1;
	fopts.by = 2;
	fopts.bz = 3;
	fopts.compute_error = GENERATE(false,true);
	fopts.subcycle = GENERATE(false,true);

	fclaw2d_domain_t *domain = create_test_domain(MPI_COMM_WORLD,&fopts);
	fclaw2d_global_store_domain(glob, domain);
	fclaw2d_options_store(glob, &fopts);

	fclaw3dx_clawpatch_options_t opts;
	memset(&opts, 0, sizeof(opts));
	opts.mx   = GENERATE(4,5,6);
	opts.my   = GENERATE(4,5,6);
	opts.mz   = GENERATE(4,5,6);
	opts.mbc  = GENERATE(1,2);
	opts.meqn = GENERATE(1,2);
	opts.maux = GENERATE(0,2);
	opts.rhs_fields = GENERATE(0,2);
	fclaw3dx_clawpatch_options_store(glob, &opts);

	fclaw2d_domain_data_new(glob->domain);
	fclaw2d_build_mode_t build_mode = GENERATE(FCLAW2D_BUILD_FOR_GHOST_AREA_COMPUTED, FCLAW2D_BUILD_FOR_UPDATE);
	CHECK(domain->blocks[0].patches[0].user == nullptr);
	fclaw2d_patch_build(glob, &domain->blocks[0].patches[0], 0, 0, &build_mode);
	CHECK(domain->blocks[0].patches[0].user != nullptr);
	fclaw3dx_clawpatch_t* cp = fclaw3dx_clawpatch_get_clawpatch(&domain->blocks[0].patches[0]);

	CHECK(cp->meqn == opts.meqn);
	CHECK(cp->mx == opts.mx);
	CHECK(cp->my == opts.my);
	CHECK(cp->mz == opts.mz);
	CHECK(cp->mbc == opts.mbc);
	CHECK(cp->manifold == fopts.manifold);
	CHECK(cp->mp != nullptr);
	CHECK(cp->registers != nullptr);

	CHECK(cp->xlower == fopts.ax);
	CHECK(cp->ylower == fopts.ay);
	CHECK(cp->zlower == fopts.az);
	CHECK(cp->xupper == fopts.bx);
	CHECK(cp->yupper == fopts.by);
	CHECK(cp->zupper == fopts.bz);
	CHECK(cp->dx == Approx((cp->xupper-cp->xlower)/opts.mx));
	CHECK(cp->dy == Approx((cp->yupper-cp->ylower)/opts.my));
	CHECK(cp->dz == Approx((cp->zupper-cp->zlower)/opts.mz));

	//BOX DIEMSIONS

	CHECK_BOX_DIMENSIONS(cp->griddata, opts.mbc, opts.mx, opts.my, opts.mz, opts.meqn);
	if(build_mode == FCLAW2D_BUILD_FOR_UPDATE){
		CHECK_BOX_DIMENSIONS(cp->griddata_last, opts.mbc, opts.mx, opts.my, opts.mz, opts.meqn);
		CHECK_BOX_DIMENSIONS(cp->griddata_save, opts.mbc, opts.mx, opts.my, opts.mz, opts.meqn);
	}else{
		CHECK_BOX_EMPTY(cp->griddata_last);
		CHECK_BOX_EMPTY(cp->griddata_save);
	}
	if(fopts.subcycle){
	    CHECK_BOX_DIMENSIONS(cp->griddata_time_interpolated, opts.mbc, opts.mx, opts.my, opts.mz, opts.meqn);
	}else{
		CHECK_BOX_EMPTY(cp->griddata_time_interpolated);
	}
	if(fopts.compute_error) {
		CHECK_BOX_DIMENSIONS(cp->griderror, opts.mbc, opts.mx, opts.my, opts.mz, opts.meqn);
		CHECK_BOX_DIMENSIONS(cp->exactsolution, opts.mbc, opts.mx, opts.my, opts.mz, opts.meqn);
	}else{
		CHECK_BOX_EMPTY(cp->griderror);
		CHECK_BOX_EMPTY(cp->exactsolution);
	}
	if(opts.rhs_fields == 0){
		CHECK_BOX_EMPTY(cp->rhs);
	}else{
		CHECK_BOX_DIMENSIONS(cp->rhs, opts.mbc, opts.mx, opts.my, opts.mz, opts.rhs_fields);
	}
	if(opts.rhs_fields == 0 || !fopts.compute_error) {
		CHECK_BOX_EMPTY(cp->elliptic_error);
		CHECK_BOX_EMPTY(cp->elliptic_soln);
	}else{
		CHECK_BOX_DIMENSIONS(cp->elliptic_error, opts.mbc, opts.mx, opts.my, opts.mz, opts.rhs_fields);
		CHECK_BOX_DIMENSIONS(cp->elliptic_soln, opts.mbc, opts.mx, opts.my, opts.mz, opts.rhs_fields);
	}

	fclaw2d_patch_data_delete(glob, &domain->blocks[0].patches[0]);
	fclaw2d_global_destroy(glob);
}