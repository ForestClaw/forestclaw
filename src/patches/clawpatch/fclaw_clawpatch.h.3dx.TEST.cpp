/*
Copyright (c) 2012-2022 Carsten Burstedde, Donna Calhoun, Scott Aiton
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

#include <fclaw_clawpatch.hpp>
#include <fclaw_clawpatch.h>
#include <fclaw_clawpatch_options.h>
#include <fclaw3dx_clawpatch46_fort.h>
#include <fclaw_clawpatch_output_ascii.h>
#include <fclaw_global.h>
#include <fclaw_domain.h>
#include <fclaw_patch.h>
#include <fclaw_convenience.h>
#include <fclaw2d_map.h>
#include <fclaw3d_metric.hpp>
#include <fclaw3d_metric.h>
#include <fclaw_options.h>
#include <test.hpp>
#include <test/test.hpp>
#include <fstream>
#include <bitset>

#include <fclaw_forestclaw.h>

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

namespace{
struct SinglePatchDomain {
    fclaw_global_t* glob;
    fclaw_options_t fopts;
    fclaw_domain_t *domain;
    fclaw2d_map_context_t* map;
    fclaw_clawpatch_options_t* opts;

    SinglePatchDomain(){
        glob = fclaw_global_new();
        opts = fclaw_clawpatch_options_new(3);
        memset(&fopts, 0, sizeof(fopts));
        fopts.mi=1;
        fopts.mj=1;
        fopts.minlevel=0;
        fopts.manifold=false;
        fopts.bx = 1;
        fopts.by = 2;
        fopts.bz = 3;
        fopts.compute_error = true;
        fopts.subcycle = true;
        fclaw_options_store(glob, &fopts);

        opts->d3->mx   = 5;
        opts->d3->my   = 6;
        opts->d3->mz   = 7;
        opts->mbc  = 2;
        opts->meqn = 1;
        opts->maux = 1;
        opts->rhs_fields = 1;
        fclaw_clawpatch_options_store(glob, opts);

        domain = fclaw_domain_new_unitsquare(sc_MPI_COMM_WORLD, 0);
        fclaw_global_store_domain(glob, domain);

        map = fclaw2d_map_new_nomap();
        fclaw2d_map_store(glob, map);

        fclaw_vtables_initialize(glob);
        fclaw_clawpatch_vtable_initialize(glob, 4);

        fclaw_domain_data_new(glob->domain);
    }
    void setup(){
        fclaw_build_mode_t build_mode = FCLAW_BUILD_FOR_UPDATE;
        fclaw_patch_build(glob, &domain->blocks[0].patches[0], 0, 0, &build_mode);
    }
    ~SinglePatchDomain(){
        fclaw_patch_data_delete(glob, &domain->blocks[0].patches[0]);
        fclaw_clawpatch_options_destroy(opts);
        fclaw_domain_destroy(domain);
        fclaw2d_map_destroy(map);
        fclaw_global_destroy(glob);
    }
};
struct QuadDomain {
    fclaw_global_t* glob;
    fclaw_options_t fopts;
    fclaw_domain_t *domain;
    fclaw2d_map_context_t* map;
    fclaw_clawpatch_options_t* opts;

    QuadDomain(){
        glob = fclaw_global_new();
        opts = fclaw_clawpatch_options_new(3);
        memset(&fopts, 0, sizeof(fopts));
        fopts.mi=1;
        fopts.mj=1;
        fopts.minlevel=1;
        fopts.maxlevel=1;
        fopts.manifold=false;
        fopts.bx = 1;
        fopts.by = 2;
        fopts.bz = 3;
        fopts.compute_error = true;
        fopts.subcycle = true;
        fclaw_options_store(glob, &fopts);

        opts->d3->mx   = 5;
        opts->d3->my   = 6;
        opts->d3->mz   = 7;
        opts->mbc  = 2;
        opts->meqn = 1;
        opts->maux = 1;
        opts->rhs_fields = 1;
        fclaw_clawpatch_options_store(glob, opts);

        domain = fclaw_domain_new_unitsquare(sc_MPI_COMM_WORLD, 1);
        fclaw_global_store_domain(glob, domain);

        map = fclaw2d_map_new_nomap();
        fclaw2d_map_store(glob, map);

        fclaw_vtables_initialize(glob);
        fclaw_clawpatch_vtable_initialize(glob, 4);

        fclaw_domain_data_new(glob->domain);
    }
    void setup(){
        fclaw_build_mode_t build_mode = FCLAW_BUILD_FOR_UPDATE;
        fclaw_patch_build(glob, &domain->blocks[0].patches[0], 0, 0, &build_mode);
        fclaw_patch_build(glob, &domain->blocks[0].patches[1], 0, 1, &build_mode);
        fclaw_patch_build(glob, &domain->blocks[0].patches[2], 0, 2, &build_mode);
        fclaw_patch_build(glob, &domain->blocks[0].patches[3], 0, 3, &build_mode);
    }
    ~QuadDomain(){
        fclaw_patch_data_delete(glob, &domain->blocks[0].patches[0]);
        fclaw_patch_data_delete(glob, &domain->blocks[0].patches[1]);
        fclaw_patch_data_delete(glob, &domain->blocks[0].patches[2]);
        fclaw_patch_data_delete(glob, &domain->blocks[0].patches[3]);
        fclaw_clawpatch_options_destroy(opts);
        fclaw_domain_destroy(domain);
        fclaw2d_map_destroy(map);
        fclaw_global_destroy(glob);
    }
};
}
TEST_CASE("3dx fclaw_clawpatch_vtable_initialize")
{
	fclaw_domain_t* domain = fclaw_domain_new_unitsquare(sc_MPI_COMM_WORLD, 1);
    fclaw_global_t* glob = fclaw_global_new();
	fclaw_global_store_domain(glob, domain);

    fclaw_clawpatch_options_t* opts = fclaw_clawpatch_options_new(3);
    fclaw_clawpatch_options_store(glob, opts);

    fclaw_vtables_initialize(glob);

    fclaw_clawpatch_vtable_initialize(glob, 4);

    fclaw_clawpatch_vtable_t * clawpatch_vt = fclaw_clawpatch_vt(glob);

    CHECK(clawpatch_vt->set_user_data             == NULL);

    //ghost filling
    CHECK(clawpatch_vt->d3->fort_copy_face              == &FCLAW3DX_CLAWPATCH46_FORT_COPY_FACE);
    CHECK(clawpatch_vt->d3->fort_average_face           == &FCLAW3DX_CLAWPATCH46_FORT_AVERAGE_FACE);
    CHECK(clawpatch_vt->d3->fort_interpolate_face       == &FCLAW3DX_CLAWPATCH46_FORT_INTERPOLATE_FACE);
    CHECK(clawpatch_vt->d3->fort_copy_corner            == &FCLAW3DX_CLAWPATCH46_FORT_COPY_CORNER);
    CHECK(clawpatch_vt->d3->fort_average_corner         == &FCLAW3DX_CLAWPATCH46_FORT_AVERAGE_CORNER);
    CHECK(clawpatch_vt->d3->fort_interpolate_corner     == &FCLAW3DX_CLAWPATCH46_FORT_INTERPOLATE_CORNER);

    //regridding
    CHECK(clawpatch_vt->d3->fort_tag4refinement         == &FCLAW3D_CLAWPATCH46_FORT_TAG4REFINEMENT);
    CHECK(clawpatch_vt->d3->fort_tag4coarsening_3dx     == &FCLAW3DX_CLAWPATCH46_FORT_TAG4COARSENING);
    CHECK(clawpatch_vt->d3->fort_tag4coarsening         == NULL);
    CHECK(clawpatch_vt->d3->fort_user_exceeds_threshold == NULL);
    CHECK(clawpatch_vt->d3->fort_interpolate2fine       == &FCLAW3DX_CLAWPATCH46_FORT_INTERPOLATE2FINE);
    CHECK(clawpatch_vt->d3->fort_average2coarse         == &FCLAW3DX_CLAWPATCH46_FORT_AVERAGE2COARSE);

    //ascii output
    CHECK(clawpatch_vt->time_header_ascii               == &fclaw_clawpatch_time_header_ascii);
    CHECK(clawpatch_vt->fort_header_ascii               == &FCLAW3D_CLAWPATCH46_FORT_HEADER_ASCII);
    CHECK(clawpatch_vt->cb_output_ascii                 == &cb_clawpatch_output_ascii);
    CHECK(clawpatch_vt->d3->fort_output_ascii           == &FCLAW3D_CLAWPATCH46_FORT_OUTPUT_ASCII);

    //time interpolation
    CHECK(clawpatch_vt->d3->fort_timeinterp             == &FCLAW3D_CLAWPATCH46_FORT_TIMEINTERP);

    //ghot packing
    CHECK(clawpatch_vt->d3->fort_local_ghost_pack       == &FCLAW3D_CLAWPATCH46_FORT_LOCAL_GHOST_PACK);
    CHECK(clawpatch_vt->local_ghost_pack_aux            == NULL);

    //diagnostics
    CHECK(clawpatch_vt->conservation_check              == NULL);
    CHECK(clawpatch_vt->compute_error                   == NULL);
    CHECK(clawpatch_vt->d3->fort_compute_patch_error    == NULL);
    CHECK(clawpatch_vt->d3->fort_conservation_check     == NULL);
    CHECK(clawpatch_vt->d3->fort_compute_error_norm     == NULL);
    CHECK(clawpatch_vt->d3->fort_compute_patch_area     == NULL);

    CHECK(clawpatch_vt->is_set                      == 1);

    fclaw_patch_vtable_t * patch_vt = fclaw_patch_vt(glob);
    //create delete build
    //TODO document patch_vt and expose these as part to public api
    CHECK(patch_vt->patch_new                      != NULL);
    CHECK(patch_vt->patch_delete                   != NULL);
    CHECK(patch_vt->build                          != NULL);
    CHECK(patch_vt->build_from_fine                != NULL);
    CHECK(patch_vt->setup                          == NULL);

    fclaw_domain_destroy(domain);
    fclaw_global_destroy(glob);
}
TEST_CASE("3dx fclaw_clawpatch patch_build")
{
    for(const int& compute_error : {0,1})
    for(const int& subcycle : {0,1})
    for(const int& mx   : {4,5,6})
    for(const int& my   : {4,5,6})
    for(const int& mz   : {4,5,6})
    for(const int& mbc  : {1,2})
    for(const int& meqn : {1,2})
    for(const int& maux : {0,2})
    for(const int& rhs_fields : {0,2})
    for(fclaw_build_mode_t build_mode : {FCLAW_BUILD_FOR_GHOST_AREA_COMPUTED, FCLAW_BUILD_FOR_UPDATE})
    {
        fclaw_global_t* glob = fclaw_global_new();

        fclaw_options_t fopts;
        memset(&fopts, 0, sizeof(fopts));
        fopts.mi=1;
        fopts.mj=1;
        fopts.minlevel=0;
        fopts.manifold=false;
        fopts.bx = 1;
        fopts.by = 2;
        fopts.bz = 3;
        fopts.compute_error = compute_error;
        fopts.subcycle = subcycle;
        fclaw_options_store(glob, &fopts);

        fclaw_clawpatch_options_t* opts = fclaw_clawpatch_options_new(3);
        opts->d3->mx     = mx;
        opts->d3->my     = my;
        opts->d3->mz     = mz;
        opts->mbc        = mbc;
        opts->meqn       = meqn;
        opts->maux       = maux;
        opts->rhs_fields = rhs_fields;
        fclaw_clawpatch_options_store(glob, opts);

	    fclaw_domain_t* domain = fclaw_domain_new_unitsquare(sc_MPI_COMM_WORLD, 0);
	    fclaw_global_store_domain(glob, domain);

        fclaw2d_map_context_t* map = fclaw2d_map_new_nomap();
        fclaw2d_map_store(glob, map);

        fclaw_vtables_initialize(glob);
        fclaw_clawpatch_vtable_initialize(glob, 4);

        fclaw_domain_data_new(glob->domain);
        CHECK(domain->blocks[0].patches[0].user == nullptr);
        fclaw_patch_build(glob, &domain->blocks[0].patches[0], 0, 0, &build_mode);
        CHECK(domain->blocks[0].patches[0].user != nullptr);

        fclaw_clawpatch_t* cp = fclaw2d_clawpatch_get_clawpatch(&domain->blocks[0].patches[0]);

        CHECK(cp->meqn == opts->meqn);
        CHECK(cp->d3->mx == opts->d3->mx);
        CHECK(cp->d3->my == opts->d3->my);
        CHECK(cp->d3->mz == opts->d3->mz);
        CHECK(cp->mbc == opts->mbc);
        CHECK(cp->manifold == fopts.manifold);
        CHECK(cp->d3->mp != nullptr);

        CHECK(cp->d3->xlower == fopts.ax);
        CHECK(cp->d3->ylower == fopts.ay);
        CHECK(cp->d3->zlower == fopts.az);
        CHECK(cp->d3->xupper == fopts.bx);
        CHECK(cp->d3->yupper == fopts.by);
        CHECK(cp->d3->zupper == fopts.bz);
        CHECK(cp->d3->dx == doctest::Approx((cp->d3->xupper-cp->d3->xlower)/opts->d3->mx));
        CHECK(cp->d3->dy == doctest::Approx((cp->d3->yupper-cp->d3->ylower)/opts->d3->my));
        CHECK(cp->d3->dz == doctest::Approx((cp->d3->zupper-cp->d3->zlower)/opts->d3->mz));

        //BOX DIEMSIONS

        CHECK_BOX_DIMENSIONS(cp->griddata, opts->mbc, opts->d3->mx, opts->d3->my, opts->d3->mz, opts->meqn);
        if(build_mode == FCLAW_BUILD_FOR_UPDATE){
            CHECK_BOX_DIMENSIONS(cp->griddata_last, opts->mbc, opts->d3->mx, opts->d3->my, opts->d3->mz, opts->meqn);
            CHECK_BOX_DIMENSIONS(cp->griddata_save, opts->mbc, opts->d3->mx, opts->d3->my, opts->d3->mz, opts->meqn);
        }else{
            CHECK_BOX_EMPTY(cp->griddata_last);
            CHECK_BOX_EMPTY(cp->griddata_save);
        }
        if(fopts.subcycle){
            CHECK_BOX_DIMENSIONS(cp->griddata_time_interpolated, opts->mbc, opts->d3->mx, opts->d3->my, opts->d3->mz, opts->meqn);
        }else{
            CHECK_BOX_EMPTY(cp->griddata_time_interpolated);
        }
        if(fopts.compute_error) {
            CHECK_BOX_DIMENSIONS(cp->griderror, opts->mbc, opts->d3->mx, opts->d3->my, opts->d3->mz, opts->meqn);
            CHECK_BOX_DIMENSIONS(cp->exactsolution, opts->mbc, opts->d3->mx, opts->d3->my, opts->d3->mz, opts->meqn);
        }else{
            CHECK_BOX_EMPTY(cp->griderror);
            CHECK_BOX_EMPTY(cp->exactsolution);
        }
        if(opts->rhs_fields == 0){
            CHECK_BOX_EMPTY(cp->rhs);
        }else{
            CHECK_BOX_DIMENSIONS(cp->rhs, opts->mbc, opts->d3->mx, opts->d3->my, opts->d3->mz, opts->rhs_fields);
        }
        if(opts->rhs_fields == 0 || !fopts.compute_error) {
            CHECK_BOX_EMPTY(cp->elliptic_error);
            CHECK_BOX_EMPTY(cp->elliptic_soln);
        }else{
            CHECK_BOX_DIMENSIONS(cp->elliptic_error, opts->mbc, opts->d3->mx, opts->d3->my, opts->d3->mz, opts->rhs_fields);
            CHECK_BOX_DIMENSIONS(cp->elliptic_soln, opts->mbc, opts->d3->mx, opts->d3->my, opts->d3->mz, opts->rhs_fields);
        }

        fclaw_patch_data_delete(glob, &domain->blocks[0].patches[0]);
        fclaw_domain_destroy(domain);
        fclaw2d_map_destroy(map);
        fclaw_global_destroy(glob);
    }
}

TEST_CASE("3dx fclaw_clawpatch save_step")
{
    SinglePatchDomain test_data;
    test_data.setup();

    fclaw_clawpatch_t* cp = fclaw2d_clawpatch_get_clawpatch(&test_data.domain->blocks[0].patches[0]);
    cp->griddata.dataPtr()[0] = 1234;
    fclaw_patch_save_step(test_data.glob,&test_data.domain->blocks[0].patches[0]);
    CHECK(cp->griddata_save.dataPtr()[0] == 1234);
}

TEST_CASE("3dx fclaw_clawpatch_save_current_step")
{
    SinglePatchDomain test_data;
    test_data.setup();

    //CHECK
    fclaw_clawpatch_t* cp = fclaw2d_clawpatch_get_clawpatch(&test_data.domain->blocks[0].patches[0]);
    cp->griddata.dataPtr()[0] = 1234;
    fclaw_clawpatch_save_current_step(test_data.glob,&test_data.domain->blocks[0].patches[0]);
    CHECK(cp->griddata_last.dataPtr()[0] == 1234);
}

TEST_CASE("3dx fclaw_clawpatch restore_step")
{
    SinglePatchDomain test_data;
    test_data.setup();

    //CHECK
    fclaw_clawpatch_t* cp = fclaw2d_clawpatch_get_clawpatch(&test_data.domain->blocks[0].patches[0]);
    cp->griddata_save.dataPtr()[0] = 1234;
    fclaw_patch_restore_step(test_data.glob,&test_data.domain->blocks[0].patches[0]);
    CHECK(cp->griddata.dataPtr()[0] == 1234);
}

#if 0
TEST_CASE("3dx fclaw_clawpatch get_metric_patch")
{
    SinglePatchDomain test_data;
    test_data.setup();

    //CHECK
    fclaw_clawpatch_t* cp = fclaw2d_clawpatch_get_clawpatch(&test_data.domain->blocks[0].patches[0]);
    CHECK(fclaw3d_patch_metric_patch(test_data.glob, &test_data.domain->blocks[0].patches[0]) == cp->mp);
}
#endif


#if 0
TEST_CASE("3dx fclaw_clawpatch_get_metric_patch")
{
    SinglePatchDomain test_data;
    test_data.setup();

    //CHECK
    fclaw_clawpatch_t* cp = 
        fclaw2d_clawpatch_get_clawpatch(&test_data.domain->blocks[0].patches[0]);
    CHECK(fclaw3d_clawpatch_get_metric_patch(&test_data.domain->blocks[0].patches[0]) == cp->mp);
}
#endif

#if 0
TEST_CASE("3dx fclaw_clawpatch_get_volume")
{
    SinglePatchDomain test_data;
    test_data.setup();

    //CHECK
    fclaw3d_metric_patch_t* mp = fclaw3dx_clawpatch_get_metric_patch(&test_data.domain->blocks[0].patches[0]);
    CHECK(fclaw3d_clawpatch_get_volume(test_data.glob, 
                &test_data.domain->blocks[0].patches[0]) == mp->volume.dataPtr());
}
#endif

TEST_CASE("3dx fclaw_clawpatch_grid_data")
{
    for(int mx   : {4,5,6})
    for(int my   : {4,5,6})
    for(int mz   : {4,5,6})
    for(int mbc  : {1,2})
    {

        SinglePatchDomain test_data;
        test_data.opts->d3->mx   = mx;
        test_data.opts->d3->my   = my;
        test_data.opts->d3->mz   = mz;
        test_data.opts->mbc  = mbc;
        test_data.setup();

        //CHECK
        fclaw_clawpatch_t* cp = fclaw2d_clawpatch_get_clawpatch(&test_data.domain->blocks[0].patches[0]);
        int mx_out,my_out,mz_out,mbc_out;
        double xlower,ylower,zlower,dx,dy,dz;
        fclaw3d_clawpatch_grid_data(test_data.glob, &test_data.domain->blocks[0].patches[0],
                                     &mx_out, &my_out, &mz_out, &mbc_out, &xlower, &ylower, &zlower, &dx, &dy, &dz);

        CHECK(mx_out == test_data.opts->d3->mx);
        CHECK(my_out == test_data.opts->d3->my);
        CHECK(mz_out == test_data.opts->d3->mz);
        CHECK(mbc_out == test_data.opts->mbc);
        CHECK(xlower == cp->d3->xlower);
        CHECK(ylower == cp->d3->ylower);
        CHECK(zlower == cp->d3->zlower);
        CHECK(dx == cp->d3->dx);
        CHECK(dy == cp->d3->dy);
        CHECK(dz == cp->d3->dz);

    }
}

TEST_CASE("3dx fclaw_clawpatch_aux_data")
{
    SinglePatchDomain test_data;
    test_data.setup();

    //CHECK
    fclaw_clawpatch_t* cp = fclaw2d_clawpatch_get_clawpatch(&test_data.domain->blocks[0].patches[0]);
    double* aux;
    int maux;
    fclaw_clawpatch_aux_data(test_data.glob, &test_data.domain->blocks[0].patches[0], &aux, &maux);

    CHECK(aux == cp->aux.dataPtr());
    CHECK(maux == test_data.opts->maux);
}

TEST_CASE("3dx fclaw_clawpatch_soln_data")
{
    SinglePatchDomain test_data;
    test_data.setup();

    //CHECK
    fclaw_clawpatch_t* cp = fclaw2d_clawpatch_get_clawpatch(&test_data.domain->blocks[0].patches[0]);
    double* q;
    int meqn;
    fclaw_clawpatch_soln_data(test_data.glob, &test_data.domain->blocks[0].patches[0], &q, &meqn);

    CHECK(q == cp->griddata.dataPtr());
    CHECK(meqn == test_data.opts->meqn);
}

TEST_CASE("3dx fclaw_clawpatch_rhs_data")
{
    SinglePatchDomain test_data;
    test_data.setup();

    //CHECK
    fclaw_clawpatch_t* cp = fclaw2d_clawpatch_get_clawpatch(&test_data.domain->blocks[0].patches[0]);
    double* rhs;
    int mfields;
    fclaw_clawpatch_rhs_data(test_data.glob, &test_data.domain->blocks[0].patches[0], &rhs, &mfields);

    CHECK(rhs == cp->rhs.dataPtr());
    CHECK(mfields == test_data.opts->rhs_fields);
}

TEST_CASE("3dx fclaw_clawpatch_elliptic_error_data")
{
    SinglePatchDomain test_data;
    test_data.setup();

    //CHECK
    fclaw_clawpatch_t* cp = fclaw2d_clawpatch_get_clawpatch(&test_data.domain->blocks[0].patches[0]);
    double* rhs;
    int mfields;
    fclaw_clawpatch_elliptic_error_data(test_data.glob, &test_data.domain->blocks[0].patches[0], &rhs, &mfields);

    CHECK(rhs == cp->elliptic_error.dataPtr());
    CHECK(mfields == test_data.opts->rhs_fields);
}

TEST_CASE("3dx fclaw_clawpatch_elliptic_soln_data")
{
    SinglePatchDomain test_data;
    test_data.setup();

    //CHECK
    fclaw_clawpatch_t* cp = fclaw2d_clawpatch_get_clawpatch(&test_data.domain->blocks[0].patches[0]);
    double* rhs;
    int mfields;
    fclaw_clawpatch_elliptic_soln_data(test_data.glob, &test_data.domain->blocks[0].patches[0], &rhs, &mfields);

    CHECK(rhs == cp->elliptic_soln.dataPtr());
    CHECK(mfields == test_data.opts->rhs_fields);
}

TEST_CASE("3dx fclaw_clawpatch_get_q")
{
    SinglePatchDomain test_data;
    test_data.setup();

    //CHECK
    fclaw_clawpatch_t* cp = fclaw2d_clawpatch_get_clawpatch(&test_data.domain->blocks[0].patches[0]);

    CHECK(fclaw_clawpatch_get_q(test_data.glob,&test_data.domain->blocks[0].patches[0]) == cp->griddata.dataPtr());
}

TEST_CASE("3dx fclaw_clawpatch_get_error")
{
    SinglePatchDomain test_data;
    test_data.setup();

    //CHECK
    fclaw_clawpatch_t* cp = fclaw2d_clawpatch_get_clawpatch(&test_data.domain->blocks[0].patches[0]);

    CHECK(fclaw_clawpatch_get_error(test_data.glob,&test_data.domain->blocks[0].patches[0]) == cp->griderror.dataPtr());
}

TEST_CASE("3dx fclaw_clawpatch_get_exact_soln")
{
    SinglePatchDomain test_data;
    test_data.setup();

    //CHECK
    fclaw_clawpatch_t* cp = fclaw2d_clawpatch_get_clawpatch(&test_data.domain->blocks[0].patches[0]);

    CHECK(fclaw_clawpatch_get_exactsoln(test_data.glob,&test_data.domain->blocks[0].patches[0]) == cp->exactsolution.dataPtr());
}
TEST_CASE("3dx fclaw_clawpatch_timesync_data")
{
    for(int time_interp : {true,false})
    {
        SinglePatchDomain test_data;
        test_data.setup();

        //CHECK
        fclaw_clawpatch_t* cp = fclaw2d_clawpatch_get_clawpatch(&test_data.domain->blocks[0].patches[0]);
        double* q;
        int meqn;
        fclaw_clawpatch_timesync_data(test_data.glob, &test_data.domain->blocks[0].patches[0], time_interp, &q, &meqn);

        if(time_interp){
            CHECK(q == cp->griddata_time_interpolated.dataPtr());
        } else {
            CHECK(q == cp->griddata.dataPtr());
        }
        CHECK(meqn == test_data.opts->meqn);
    }
}
TEST_CASE("3dx fclaw_clawpatch_get_q_timesync")
{
    for(int time_interp : {true,false})
    {
        SinglePatchDomain test_data;
        test_data.setup();

        //CHECK
        fclaw_clawpatch_t* cp = fclaw2d_clawpatch_get_clawpatch(&test_data.domain->blocks[0].patches[0]);

        if(time_interp){
            CHECK(fclaw_clawpatch_get_q_timesync(test_data.glob,&test_data.domain->blocks[0].patches[0],time_interp) == cp->griddata_time_interpolated.dataPtr());
        } else {
            CHECK(fclaw_clawpatch_get_q_timesync(test_data.glob,&test_data.domain->blocks[0].patches[0],time_interp) == cp->griddata.dataPtr());
        }
    }
}
TEST_CASE("3dx fclaw_clawpatch user_data")
{
    SinglePatchDomain test_data;
    test_data.setup();

    //CHECK
    void* user_data = (void *) 1234;
    fclaw_clawpatch_set_user_data(test_data.glob, &test_data.domain->blocks[0].patches[0],user_data);
    CHECK(fclaw_clawpatch_get_user_data(test_data.glob,&test_data.domain->blocks[0].patches[0]) == user_data);
}
TEST_CASE("3dx fclaw_clawpatch solver_data")
{
    SinglePatchDomain test_data;
    test_data.setup();

    //CHECK
    void* user_data = (void *) 1234;
    fclaw_clawpatch_set_solver_data(test_data.glob, &test_data.domain->blocks[0].patches[0],user_data);
    CHECK(fclaw_clawpatch_get_solver_data(test_data.glob,&test_data.domain->blocks[0].patches[0]) == user_data);
}

TEST_CASE("3dx fclaw_clawpatch_size")
{
    for(int mx   : {4,5,6})
    for(int my   : {4,5,6})
    for(int mz   : {4,5,6})
    for(int mbc  : {1,2})
    for(int meqn : {1,2})
    {
        SinglePatchDomain test_data;
        test_data.opts->d3->mx   = mx;
        test_data.opts->d3->my   = my;
        test_data.opts->d3->mz   = mz;
        test_data.opts->mbc  = mbc;
        test_data.opts->meqn = meqn;
        test_data.setup();

        //CHECK
        CHECK(fclaw_clawpatch_size(test_data.glob) == (size_t) (mx+2*mbc)*(my+2*mbc)*(mz+2*mbc)*meqn);
    }
}

#if 0
TEST_CASE("3dx fclaw_clawpatch_metric_scalar")
{
    SinglePatchDomain test_data;
    test_data.setup();

    double *mp_volume, *mp_faceareas;
    double *volume, *faceareas;

    /* DAC : What is the right test here? */
    fclaw3d_metric_patch_scalar(test_data.glob, 
                                &test_data.domain->blocks[0].patches[0], 
                                &mp_volume, &mp_faceareas);

    fclaw3dx_clawpatch_metric_scalar(test_data.glob, 
                                     &test_data.domain->blocks[0].patches[0], 
                                     &volume, &faceareas);

    CHECK(mp_volume == volume);
    CHECK(mp_faceareas == faceareas);
}
#endif

#if 0
/* DAC : Fix this test for fclaw3d_metric? */
TEST_CASE("3dx fclaw_clawpatch_metric_vector")
{
    fclaw_global_t* glob = fclaw_global_new(); 
    fclaw_vtables_initialize(glob);

    SinglePatchDomain test_data;
    test_data.setup();

    double *mp_xnormals, *mp_ynormals, *mp_xtangents, *mp_ytangents, *mp_curvature;
    double *xnormals, *ynormals, *xtangents, *ytangents, *curvature;

    fclaw3d_metric_patch_vector(test_data.glob, &test_data.domain->blocks[0].patches[0], 
                                &mp_xnormals, &mp_ynormals,
                                &mp_xtangents, &mp_ytangents,
                                &mp_curvature);

    fclaw3dx_clawpatch_metric_vector(test_data.glob, &test_data.domain->blocks[0].patches[0], 
                                     &xnormals, &ynormals,
                                      &xtangents, &ytangents,
                                      &curvature);

    CHECK(mp_xnormals == xnormals);
    CHECK(mp_ynormals == ynormals);
    CHECK(mp_xtangents == xtangents);
    CHECK(mp_ytangents == ytangents);
    CHECK(mp_curvature == curvature);

    fclaw_global_destroy(glob);
}
#endif

#if 0
TEST_CASE("3dx fclaw_clawpatch_metric_data")
{
    SinglePatchDomain test_data;
    test_data.setup();

    double *mp_xp, *mp_yp, *mp_zp, *mp_xd, *mp_yd, *mp_zd, *mp_area;
    double *xp, *yp, *zp, *xd, *yd, *zd, *area;

    fclaw3d_metric_patch_mesh_data(test_data.glob, &test_data.domain->blocks[0].patches[0], 
                                   &mp_xp, &mp_yp, &mp_zp,
                                   &mp_xd, &mp_yd, &mp_zd,
                                   &mp_area);

    fclaw3dx_clawpatch_metric_data(test_data.glob, &test_data.domain->blocks[0].patches[0], 
                                   &xp, &yp, &zp,
                                   &xd, &yd, &zd,
                                   &area);

    CHECK(mp_xp == xp);
    CHECK(mp_yp == yp);
    CHECK(mp_zp == zp);
    CHECK(mp_xd == xd);
    CHECK(mp_yd == yd);
    CHECK(mp_zd == zd);
    CHECK(mp_area == area);
}
#endif

#if 0
TEST_CASE("3dx fclaw_clawpatch_metric_data2")
{
    SinglePatchDomain test_data;
    test_data.setup();

    double *mp_xnormals, *mp_ynormals, *mp_xtangents, *mp_ytangents, *mp_surfnormals, *mp_edgelengths, *mp_curvature;
    double *xnormals, *ynormals, *xtangents, *ytangents, *surfnormals, *edgelengths, *curvature;

    fclaw3d_metric_patch_mesh_data2(test_data.glob, &test_data.domain->blocks[0].patches[0], 
                                    &mp_xnormals, &mp_ynormals,
                                    &mp_xtangents, &mp_ytangents,
                                    &mp_surfnormals, &mp_edgelengths, &mp_curvature);

    fclaw3dx_clawpatch_metric_data2(test_data.glob, &test_data.domain->blocks[0].patches[0], 
                                     &xnormals, &ynormals,
                                     &xtangents, &ytangents,
                                     &surfnormals, &edgelengths, &curvature);

    CHECK(mp_xnormals == xnormals);
    CHECK(mp_ynormals == ynormals);
    CHECK(mp_xtangents == xtangents);
    CHECK(mp_ytangents == ytangents);
    CHECK(mp_surfnormals == surfnormals);
    CHECK(mp_edgelengths == edgelengths);
    CHECK(mp_curvature == curvature);
}
#endif

namespace{
    double timeinterp_alpha;
    int timeinterp_mint;
    fclaw_clawpatch_t* timeinterp_cp;
}


TEST_CASE("3dx fclaw_clawpatch setup_timeinterp")
{
    SinglePatchDomain test_data;
    test_data.opts->interp_stencil_width=2;
    timeinterp_mint = 2;
    test_data.setup();


    timeinterp_cp = fclaw2d_clawpatch_get_clawpatch(&test_data.domain->blocks[0].patches[0]);
    timeinterp_alpha = 0.90210;


    fclaw_clawpatch_vtable_t * clawpatch_vt = fclaw_clawpatch_vt(test_data.glob);
    clawpatch_vt->d3->fort_timeinterp = [] (const int *mx, const int *my, const int *mz, 
                                           const int *mbc, const int *meqn, const int *psize, 
                                           double qcurr[], double qlast[], double qinterp[], 
                                           const double *alpha, const int *ierror)
    {
        CHECK(*mx == timeinterp_cp->d3->mx);
        CHECK(*my == timeinterp_cp->d3->my);
        CHECK(*mz == timeinterp_cp->d3->mz);
        CHECK(*mbc == timeinterp_cp->mbc);
        CHECK(*meqn == timeinterp_cp->meqn);
        CHECK(*psize == (*mx)*(*my)*(*mz)
                        - (*mx-2*timeinterp_mint)*(*my-2*timeinterp_mint)*(*mz));
        CHECK(qcurr == timeinterp_cp->griddata.dataPtr());
        CHECK(qlast == timeinterp_cp->griddata_last.dataPtr());
        CHECK(qinterp == timeinterp_cp->griddata_time_interpolated.dataPtr());
        CHECK(*alpha == timeinterp_alpha);
    };

    fclaw_patch_setup_timeinterp(test_data.glob, &test_data.domain->blocks[0].patches[0], 
                                 timeinterp_alpha);
}
namespace{
    fclaw_clawpatch_t* t4r_cp;
    int t4r_tag_patch;
    int t4r_init_flag;
}
TEST_CASE("3dx fclaw_clawpatch tag4refinement")
{
    for(int tag_patch : {true,false})
    for(int init_flag : {true,false})
    {
        SinglePatchDomain test_data;
        test_data.fopts.refine_threshold = 0.90210;
        test_data.setup();


        t4r_cp = fclaw2d_clawpatch_get_clawpatch(&test_data.domain->blocks[0].patches[0]);

        fclaw_clawpatch_vtable_t * clawpatch_vt = fclaw_clawpatch_vt(test_data.glob);

        t4r_tag_patch = tag_patch;
        t4r_init_flag = init_flag;
        clawpatch_vt->d3->fort_tag4refinement = [](const int *mx, const int *my, const int *mz, 
                                                   const int *mbc, const int *meqn, 
                                                   const double *xlower, const double *ylower, const double *zlower, 
                                                   const double *dx, const double *dy, const double *dz, 
                                                   const int *blockno, 
                                                   double q[], const double *tag_threshold, const int *init_flag, 
                                                   int *tag_patch)
        {
            CHECK(*mx == t4r_cp->d3->mx);
            CHECK(*my == t4r_cp->d3->my);
            CHECK(*mz == t4r_cp->d3->mz);
            CHECK(*mbc == t4r_cp->mbc);
            CHECK(*meqn == t4r_cp->meqn);
            CHECK(*xlower == t4r_cp->d3->xlower);
            CHECK(*ylower == t4r_cp->d3->ylower);
            CHECK(*zlower == t4r_cp->d3->zlower);
            CHECK(*dx == t4r_cp->d3->dx);
            CHECK(*dy == t4r_cp->d3->dy);
            CHECK(*dz == t4r_cp->d3->dz);
            CHECK(*blockno == t4r_cp->blockno);
            CHECK(q == t4r_cp->griddata.dataPtr());
            CHECK(*tag_threshold == .90210);
            CHECK(*init_flag == t4r_init_flag);
            *tag_patch = t4r_tag_patch;
        };

        CHECK(fclaw_patch_tag4refinement(test_data.glob, &test_data.domain->blocks[0].patches[0], 0, 0, t4r_init_flag) == t4r_tag_patch);
    }
}
TEST_CASE("3dx fclaw_clawpatch tag4refinement negative refine threshold")
{
    for(int init_flag : {true,false})
    {
        SinglePatchDomain test_data;
        test_data.fopts.refine_threshold = -0.90210;
        test_data.setup();

        t4r_init_flag = init_flag;

        CHECK(fclaw_patch_tag4refinement(test_data.glob, &test_data.domain->blocks[0].patches[0], 0, 0, t4r_init_flag) == true);
    }
}
namespace{
    fclaw_clawpatch_t* t4c_cp;
    fclaw_clawpatch_t* t4c_cp1;
    fclaw_clawpatch_t* t4c_cp2;
    fclaw_clawpatch_t* t4c_cp3;
    int t4c_tag_patch;
    int t4c_init_flag;
}
TEST_CASE("3dx fclaw_clawpatch tag4coarsening")
{
    for(int tag_patch : {true,false})
    {
        t4c_tag_patch= tag_patch;
        for(int init_flag : {true,false})
        {
            t4c_init_flag = init_flag;

            QuadDomain test_data;
            test_data.fopts.coarsen_threshold = 0.90210;
            test_data.setup();


            t4c_cp = fclaw2d_clawpatch_get_clawpatch(&test_data.domain->blocks[0].patches[0]);
            t4c_cp1 = fclaw2d_clawpatch_get_clawpatch(&test_data.domain->blocks[0].patches[1]);
            t4c_cp2 = fclaw2d_clawpatch_get_clawpatch(&test_data.domain->blocks[0].patches[2]);
            t4c_cp3 = fclaw2d_clawpatch_get_clawpatch(&test_data.domain->blocks[0].patches[3]);

            fclaw_clawpatch_vtable_t * clawpatch_vt = fclaw_clawpatch_vt(test_data.glob);

            clawpatch_vt->d3->fort_tag4coarsening_3dx = [](const int *mx, const int *my, const int *mz, 
                                                       const int *mbc, const int *meqn, 
                                                       double xlower[], double ylower[], double zlower[], 
                                                       const double *dx, const double *dy, const double *dz, 
                                                       const int *blockno, 
                                                       double q0[], double q1[], double q2[], double q3[], 
                                                       const double *tag_threshold, const int *init_flag, int *tag_patch)
            {
                CHECK(*mx == t4c_cp->d3->mx);
                CHECK(*my == t4c_cp->d3->my);
                CHECK(*mz == t4c_cp->d3->mz);
                CHECK(*mbc == t4c_cp->mbc);
                CHECK(*meqn == t4c_cp->meqn);
                CHECK(*xlower == t4c_cp->d3->xlower);
                CHECK(*ylower == t4c_cp->d3->ylower);
                CHECK(*zlower == t4c_cp->d3->zlower);
                CHECK(*dx == t4c_cp->d3->dx);
                CHECK(*dy == t4c_cp->d3->dy);
                CHECK(*dz == t4c_cp->d3->dz);
                CHECK(*blockno == t4c_cp->blockno);
                CHECK(q0 == t4c_cp->griddata.dataPtr());
                CHECK(q1 == t4c_cp1->griddata.dataPtr());
                CHECK(q2 == t4c_cp2->griddata.dataPtr());
                CHECK(q3 == t4c_cp3->griddata.dataPtr());
                CHECK(*tag_threshold == .90210);
                CHECK(*init_flag == t4c_init_flag);
                *tag_patch = t4c_tag_patch;
            };

            CHECK(fclaw_patch_tag4coarsening(test_data.glob, &test_data.domain->blocks[0].patches[0], 0, 0, t4c_init_flag) == t4c_tag_patch);
        }
    }
}
TEST_CASE("3dx fclaw_clawpatch tag4coarsening negative coarsen threshold")
{
    for(int init_flag : {true,false})
    {
        t4c_init_flag = init_flag;
    
        QuadDomain test_data;
        test_data.fopts.coarsen_threshold = -0.90210;
        test_data.setup();

        CHECK(fclaw_patch_tag4coarsening(test_data.glob, &test_data.domain->blocks[0].patches[0], 0, 0, t4c_init_flag) == false);
    }
}
namespace{
    fclaw_clawpatch_t* i2f_ccp;
    fclaw_clawpatch_t* i2f_cp;
    fclaw_clawpatch_t* i2f_cp1;
    fclaw_clawpatch_t* i2f_cp2;
    fclaw_clawpatch_t* i2f_cp3;
    std::bitset<4> i2f_igrids;
    int i2f_manifold;
}
TEST_CASE("3dx fclaw_clawpatch interpolate2fine")
{
    i2f_manifold = false;

    QuadDomain fine_test_data;
    fine_test_data.fopts.manifold = i2f_manifold;
    fine_test_data.setup();
    SinglePatchDomain coarse_test_data;
    coarse_test_data.fopts.manifold = i2f_manifold;
    coarse_test_data.setup();


    i2f_ccp = fclaw2d_clawpatch_get_clawpatch(&coarse_test_data.domain->blocks[0].patches[0]);
    i2f_cp = fclaw2d_clawpatch_get_clawpatch(&fine_test_data.domain->blocks[0].patches[0]);
    i2f_cp1 = fclaw2d_clawpatch_get_clawpatch(&fine_test_data.domain->blocks[0].patches[1]);
    i2f_cp2 = fclaw2d_clawpatch_get_clawpatch(&fine_test_data.domain->blocks[0].patches[2]);
    i2f_cp3 = fclaw2d_clawpatch_get_clawpatch(&fine_test_data.domain->blocks[0].patches[3]);

    fclaw_clawpatch_vtable_t * clawpatch_vt = fclaw_clawpatch_vt(coarse_test_data.glob);

    clawpatch_vt->d3->fort_interpolate2fine = [](const int *mx, const int *my, const int *mz, 
                                                 const int *mbc, const int *meqn, 
                                                 double qcoarse[], double qfine[], double areacoarse[], double areafine[], 
                                                 const int *igrid, const int *manifold)
    {
        CHECK(*mx == i2f_cp->d3->mx);
        CHECK(*my == i2f_cp->d3->my);
        CHECK(*mz == i2f_cp->d3->mz);
        CHECK(*mbc == i2f_cp->mbc);
        CHECK(*meqn == i2f_cp->meqn);

        CHECK(qcoarse == i2f_ccp->griddata.dataPtr());
#if PATCH_DIM == 2
        CHECK(areacoarse == i2f_ccp->mp->area.dataPtr());

        if(*igrid==0){
            CHECK(qfine == i2f_cp->griddata.dataPtr());
            CHECK(areafine == i2f_cp->mp->area.dataPtr());
        }else if(*igrid==1){
            CHECK(qfine == i2f_cp1->griddata.dataPtr());
            CHECK(areafine == i2f_cp1->mp->area.dataPtr());
        }else if(*igrid==2){
            CHECK(qfine == i2f_cp2->griddata.dataPtr());
            CHECK(areafine == i2f_cp2->mp->area.dataPtr());
        }else if(*igrid==3){
            CHECK(qfine == i2f_cp3->griddata.dataPtr());
            CHECK(areafine == i2f_cp3->mp->area.dataPtr());
        }else{
            REQUIRE(false);
        }
#endif

        CHECK(i2f_igrids[*igrid] == false);
        i2f_igrids[*igrid] = true;
        CHECK(*manifold == i2f_manifold);
    };

    fclaw_patch_interpolate2fine(coarse_test_data.glob,
                                   &coarse_test_data.domain->blocks[0].patches[0],
                                   &fine_test_data.domain->blocks[0].patches[0],
                                   0, 0, 0);

    CHECK(i2f_igrids.all());
}