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

#include <fclaw3dx_clawpatch.h>
#include <fclaw3dx_clawpatch.hpp>
#include <fclaw3dx_clawpatch_options.h>
#include <fclaw3dx_clawpatch46_fort.h>
#include <fclaw3dx_clawpatch_output_ascii.h>
#include <fclaw2d_global.h>
#include <fclaw2d_domain.h>
#include <fclaw2d_patch.h>
#include <fclaw2d_convenience.h>
#include <fclaw2d_metric.hpp>
#include <fclaw2d_metric.h>
#include <fclaw2d_options.h>
#include <test/doctest.h>
#include <test/test.hpp>
#include <fstream>
#include <bitset>

#include <fclaw2d_forestclaw.h>

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
    fclaw2d_global_t* glob;
    fclaw_options_t fopts;
    fclaw2d_domain_t *domain;
    fclaw3dx_clawpatch_options_t opts;

    SinglePatchDomain(){
        glob = fclaw2d_global_new();
        fclaw2d_vtables_initialize(glob);
        fclaw3dx_clawpatch_vtable_initialize(glob, 4);
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

        domain = create_test_domain(sc_MPI_COMM_WORLD,&fopts);
        fclaw2d_global_store_domain(glob, domain);
        fclaw2d_options_store(glob, &fopts);

        memset(&opts, 0, sizeof(opts));
        opts.mx   = 5;
        opts.my   = 6;
        opts.mz   = 7;
        opts.mbc  = 2;
        opts.meqn = 1;
        opts.maux = 1;
        opts.rhs_fields = 1;
        fclaw3dx_clawpatch_options_store(glob, &opts);

        fclaw2d_domain_data_new(glob->domain);
    }
    void setup(){
        fclaw2d_build_mode_t build_mode = FCLAW2D_BUILD_FOR_UPDATE;
        fclaw2d_patch_build(glob, &domain->blocks[0].patches[0], 0, 0, &build_mode);
    }
    ~SinglePatchDomain(){
        fclaw2d_patch_data_delete(glob, &domain->blocks[0].patches[0]);
        fclaw2d_global_destroy(glob);
    }
};
struct QuadDomain {
    fclaw2d_global_t* glob;
    fclaw_options_t fopts;
    fclaw2d_domain_t *domain;
    fclaw3dx_clawpatch_options_t opts;

    QuadDomain(){
        glob = fclaw2d_global_new();
        fclaw2d_vtables_initialize(glob);
        fclaw3dx_clawpatch_vtable_initialize(glob, 4);
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

        domain = create_test_domain(sc_MPI_COMM_WORLD,&fopts);
        fclaw2d_global_store_domain(glob, domain);
        fclaw2d_options_store(glob, &fopts);

        memset(&opts, 0, sizeof(opts));
        opts.mx   = 5;
        opts.my   = 6;
        opts.mz   = 7;
        opts.mbc  = 2;
        opts.meqn = 1;
        opts.maux = 1;
        opts.rhs_fields = 1;
        fclaw3dx_clawpatch_options_store(glob, &opts);

        fclaw2d_domain_data_new(glob->domain);
    }
    void setup(){
        fclaw2d_build_mode_t build_mode = FCLAW2D_BUILD_FOR_UPDATE;
        fclaw2d_patch_build(glob, &domain->blocks[0].patches[0], 0, 0, &build_mode);
        fclaw2d_patch_build(glob, &domain->blocks[0].patches[1], 0, 1, &build_mode);
        fclaw2d_patch_build(glob, &domain->blocks[0].patches[2], 0, 2, &build_mode);
        fclaw2d_patch_build(glob, &domain->blocks[0].patches[3], 0, 3, &build_mode);
    }
    ~QuadDomain(){
        fclaw2d_patch_data_delete(glob, &domain->blocks[0].patches[0]);
        fclaw2d_patch_data_delete(glob, &domain->blocks[0].patches[1]);
        fclaw2d_patch_data_delete(glob, &domain->blocks[0].patches[2]);
        fclaw2d_patch_data_delete(glob, &domain->blocks[0].patches[3]);
        fclaw2d_global_destroy(glob);
    }
};
}
TEST_CASE("fclaw3dx_clawpatch_vtable_initialize")
{
    fclaw2d_global_t* glob = fclaw2d_global_new();
    fclaw2d_vtables_initialize(glob);

    fclaw3dx_clawpatch_vtable_initialize(glob, 4);

    fclaw3dx_clawpatch_vtable_t * clawpatch_vt = fclaw3dx_clawpatch_vt(glob);

    CHECK(clawpatch_vt->set_user_data             != NULL);

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
    CHECK(clawpatch_vt->conservation_check          != NULL);
    CHECK(clawpatch_vt->compute_error               != NULL);
    CHECK(clawpatch_vt->fort_compute_patch_error    != NULL);
    CHECK(clawpatch_vt->fort_conservation_check     != NULL);
    CHECK(clawpatch_vt->fort_compute_error_norm     != NULL);
    CHECK(clawpatch_vt->fort_compute_patch_area     != NULL);

    CHECK(clawpatch_vt->is_set                      == 1);

    fclaw2d_patch_vtable_t * patch_vt = fclaw2d_patch_vt(glob);
    //create delete build
    //TODO document patch_vt and expose these as part to public api
    CHECK(patch_vt->patch_new                      != NULL);
    CHECK(patch_vt->patch_delete                   != NULL);
    CHECK(patch_vt->build                          != NULL);
    CHECK(patch_vt->build_from_fine                != NULL);
    CHECK(patch_vt->setup                          == NULL);
    fclaw2d_global_destroy(glob);
}
TEST_CASE("fclaw3dx_clawpatch patch_build")
{
    for(const int& compute_error : {false,true})
    for(const int& subcycle : {false,true})
    for(const int& mx   : {4,5,6})
    for(const int& my   : {4,5,6})
    for(const int& mz   : {4,5,6})
    for(const int& mbc  : {1,2})
    for(const int& meqn : {1,2})
    for(const int& maux : {0,2})
    for(const int& rhs_fields : {0,2})
    for(fclaw2d_build_mode_t build_mode : {FCLAW2D_BUILD_FOR_GHOST_AREA_COMPUTED, FCLAW2D_BUILD_FOR_UPDATE})
    {
        fclaw2d_global_t* glob = fclaw2d_global_new(); 
        fclaw2d_vtables_initialize(glob);

        fclaw3dx_clawpatch_vtable_initialize(glob, 4);

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

        fclaw2d_domain_t *domain = create_test_domain(sc_MPI_COMM_WORLD,&fopts);
        fclaw2d_global_store_domain(glob, domain);
        fclaw2d_options_store(glob, &fopts);

        fclaw3dx_clawpatch_options_t opts;
        memset(&opts, 0, sizeof(opts));
        opts.mx   = mx;
        opts.my   = my;
        opts.mz   = mz;
        opts.mbc  = mbc;
        opts.meqn = meqn;
        opts.maux = maux;
        opts.rhs_fields = rhs_fields;
        fclaw3dx_clawpatch_options_store(glob, &opts);

        fclaw2d_domain_data_new(glob->domain);
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
        CHECK(cp->dx == doctest::Approx((cp->xupper-cp->xlower)/opts.mx));
        CHECK(cp->dy == doctest::Approx((cp->yupper-cp->ylower)/opts.my));
        CHECK(cp->dz == doctest::Approx((cp->zupper-cp->zlower)/opts.mz));

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
}

TEST_CASE("fclaw3dx_clawpatch save_step")
{
    SinglePatchDomain test_data;
    test_data.setup();

    fclaw3dx_clawpatch_t* cp = fclaw3dx_clawpatch_get_clawpatch(&test_data.domain->blocks[0].patches[0]);
    cp->griddata.dataPtr()[0] = 1234;
    fclaw2d_patch_save_step(test_data.glob,&test_data.domain->blocks[0].patches[0]);
    CHECK(cp->griddata_save.dataPtr()[0] == 1234);
}

TEST_CASE("fclaw3dx_clawpatch_save_current_step")
{
    SinglePatchDomain test_data;
    test_data.setup();

    //CHECK
    fclaw3dx_clawpatch_t* cp = fclaw3dx_clawpatch_get_clawpatch(&test_data.domain->blocks[0].patches[0]);
    cp->griddata.dataPtr()[0] = 1234;
    fclaw3dx_clawpatch_save_current_step(test_data.glob,&test_data.domain->blocks[0].patches[0]);
    CHECK(cp->griddata_last.dataPtr()[0] == 1234);
}

TEST_CASE("fclaw3dx_clawpatch restore_step")
{
    SinglePatchDomain test_data;
    test_data.setup();

    //CHECK
    fclaw3dx_clawpatch_t* cp = fclaw3dx_clawpatch_get_clawpatch(&test_data.domain->blocks[0].patches[0]);
    cp->griddata_save.dataPtr()[0] = 1234;
    fclaw2d_patch_restore_step(test_data.glob,&test_data.domain->blocks[0].patches[0]);
    CHECK(cp->griddata.dataPtr()[0] == 1234);
}

TEST_CASE("fclaw3dx_clawpatch get_metric_patch")
{
    SinglePatchDomain test_data;
    test_data.setup();

    //CHECK
    fclaw3dx_clawpatch_t* cp = fclaw3dx_clawpatch_get_clawpatch(&test_data.domain->blocks[0].patches[0]);
    CHECK(fclaw2d_patch_metric_patch(test_data.glob, &test_data.domain->blocks[0].patches[0]) == cp->mp);
}

TEST_CASE("fclaw3dx_clawpatch_get_metric_patch")
{
    SinglePatchDomain test_data;
    test_data.setup();

    //CHECK
    fclaw3dx_clawpatch_t* cp = fclaw3dx_clawpatch_get_clawpatch(&test_data.domain->blocks[0].patches[0]);
    CHECK(fclaw3dx_clawpatch_get_metric_patch(&test_data.domain->blocks[0].patches[0]) == cp->mp);
}

TEST_CASE("fclaw3dx_clawpatch_get_area")
{
    SinglePatchDomain test_data;
    test_data.setup();

    //CHECK
    fclaw2d_metric_patch_t* mp = fclaw3dx_clawpatch_get_metric_patch(&test_data.domain->blocks[0].patches[0]);
    CHECK(fclaw3dx_clawpatch_get_area(test_data.glob, &test_data.domain->blocks[0].patches[0]) == mp->area.dataPtr());
}

TEST_CASE("fclaw3dx_clawpatch_grid_data")
{
    for(int mx   : {4,5,6})
    for(int my   : {4,5,6})
    for(int mz   : {4,5,6})
    for(int mbc  : {1,2})
    {

        SinglePatchDomain test_data;
        test_data.opts.mx   = mx;
        test_data.opts.my   = my;
        test_data.opts.mz   = mz;
        test_data.opts.mbc  = mbc;
        test_data.setup();

        //CHECK
        fclaw3dx_clawpatch_t* cp = fclaw3dx_clawpatch_get_clawpatch(&test_data.domain->blocks[0].patches[0]);
        int mx_out,my_out,mz_out,mbc_out;
        double xlower,ylower,zlower,dx,dy,dz;
        fclaw3dx_clawpatch_grid_data(test_data.glob, &test_data.domain->blocks[0].patches[0],
                                     &mx_out, &my_out, &mz_out, &mbc_out, &xlower, &ylower, &zlower, &dx, &dy, &dz);

        CHECK(mx_out == test_data.opts.mx);
        CHECK(my_out == test_data.opts.my);
        CHECK(mz_out == test_data.opts.mz);
        CHECK(mbc_out == test_data.opts.mbc);
        CHECK(xlower == cp->xlower);
        CHECK(ylower == cp->ylower);
        CHECK(zlower == cp->zlower);
        CHECK(dx == cp->dx);
        CHECK(dy == cp->dy);
        CHECK(dz == cp->dz);

    }
}

TEST_CASE("fclaw3dx_clawpatch_aux_data")
{
    SinglePatchDomain test_data;
    test_data.setup();

    //CHECK
    fclaw3dx_clawpatch_t* cp = fclaw3dx_clawpatch_get_clawpatch(&test_data.domain->blocks[0].patches[0]);
    double* aux;
    int maux;
    fclaw3dx_clawpatch_aux_data(test_data.glob, &test_data.domain->blocks[0].patches[0], &aux, &maux);

    CHECK(aux == cp->aux.dataPtr());
    CHECK(maux == test_data.opts.maux);
}

TEST_CASE("fclaw3dx_clawpatch_soln_data")
{
    SinglePatchDomain test_data;
    test_data.setup();

    //CHECK
    fclaw3dx_clawpatch_t* cp = fclaw3dx_clawpatch_get_clawpatch(&test_data.domain->blocks[0].patches[0]);
    double* q;
    int meqn;
    fclaw3dx_clawpatch_soln_data(test_data.glob, &test_data.domain->blocks[0].patches[0], &q, &meqn);

    CHECK(q == cp->griddata.dataPtr());
    CHECK(meqn == test_data.opts.meqn);
}

TEST_CASE("fclaw3dx_clawpatch_rhs_data")
{
    SinglePatchDomain test_data;
    test_data.setup();

    //CHECK
    fclaw3dx_clawpatch_t* cp = fclaw3dx_clawpatch_get_clawpatch(&test_data.domain->blocks[0].patches[0]);
    double* rhs;
    int mfields;
    fclaw3dx_clawpatch_rhs_data(test_data.glob, &test_data.domain->blocks[0].patches[0], &rhs, &mfields);

    CHECK(rhs == cp->rhs.dataPtr());
    CHECK(mfields == test_data.opts.rhs_fields);
}

TEST_CASE("fclaw3dx_clawpatch_elliptic_error_data")
{
    SinglePatchDomain test_data;
    test_data.setup();

    //CHECK
    fclaw3dx_clawpatch_t* cp = fclaw3dx_clawpatch_get_clawpatch(&test_data.domain->blocks[0].patches[0]);
    double* rhs;
    int mfields;
    fclaw3dx_clawpatch_elliptic_error_data(test_data.glob, &test_data.domain->blocks[0].patches[0], &rhs, &mfields);

    CHECK(rhs == cp->elliptic_error.dataPtr());
    CHECK(mfields == test_data.opts.rhs_fields);
}

TEST_CASE("fclaw3dx_clawpatch_elliptic_soln_data")
{
    SinglePatchDomain test_data;
    test_data.setup();

    //CHECK
    fclaw3dx_clawpatch_t* cp = fclaw3dx_clawpatch_get_clawpatch(&test_data.domain->blocks[0].patches[0]);
    double* rhs;
    int mfields;
    fclaw3dx_clawpatch_elliptic_soln_data(test_data.glob, &test_data.domain->blocks[0].patches[0], &rhs, &mfields);

    CHECK(rhs == cp->elliptic_soln.dataPtr());
    CHECK(mfields == test_data.opts.rhs_fields);
}

TEST_CASE("fclaw3dx_clawpatch_get_q")
{
    SinglePatchDomain test_data;
    test_data.setup();

    //CHECK
    fclaw3dx_clawpatch_t* cp = fclaw3dx_clawpatch_get_clawpatch(&test_data.domain->blocks[0].patches[0]);

    CHECK(fclaw3dx_clawpatch_get_q(test_data.glob,&test_data.domain->blocks[0].patches[0]) == cp->griddata.dataPtr());
}

TEST_CASE("fclaw3dx_clawpatch_get_registers")
{
    SinglePatchDomain test_data;
    test_data.setup();

    //CHECK
    fclaw3dx_clawpatch_t* cp = fclaw3dx_clawpatch_get_clawpatch(&test_data.domain->blocks[0].patches[0]);

    CHECK(fclaw3dx_clawpatch_get_registers(test_data.glob,&test_data.domain->blocks[0].patches[0]) == cp->registers);
}

TEST_CASE("fclaw3dx_clawpatch_get_error")
{
    SinglePatchDomain test_data;
    test_data.setup();

    //CHECK
    fclaw3dx_clawpatch_t* cp = fclaw3dx_clawpatch_get_clawpatch(&test_data.domain->blocks[0].patches[0]);

    CHECK(fclaw3dx_clawpatch_get_error(test_data.glob,&test_data.domain->blocks[0].patches[0]) == cp->griderror.dataPtr());
}

TEST_CASE("fclaw3dx_clawpatch_get_exact_soln")
{
    SinglePatchDomain test_data;
    test_data.setup();

    //CHECK
    fclaw3dx_clawpatch_t* cp = fclaw3dx_clawpatch_get_clawpatch(&test_data.domain->blocks[0].patches[0]);

    CHECK(fclaw3dx_clawpatch_get_exactsoln(test_data.glob,&test_data.domain->blocks[0].patches[0]) == cp->exactsolution.dataPtr());
}
TEST_CASE("fclaw3dx_clawpatch_timesync_data")
{
    for(int time_interp : {true,false})
    {
        SinglePatchDomain test_data;
        test_data.setup();

        //CHECK
        fclaw3dx_clawpatch_t* cp = fclaw3dx_clawpatch_get_clawpatch(&test_data.domain->blocks[0].patches[0]);
        double* q;
        int meqn;
        fclaw3dx_clawpatch_timesync_data(test_data.glob, &test_data.domain->blocks[0].patches[0], time_interp, &q, &meqn);

        if(time_interp){
            CHECK(q == cp->griddata_time_interpolated.dataPtr());
        } else {
            CHECK(q == cp->griddata.dataPtr());
        }
        CHECK(meqn == test_data.opts.meqn);
    }
}
TEST_CASE("fclaw3dx_clawpatch_get_q_timesync")
{
    for(int time_interp : {true,false})
    {
        SinglePatchDomain test_data;
        test_data.setup();

        //CHECK
        fclaw3dx_clawpatch_t* cp = fclaw3dx_clawpatch_get_clawpatch(&test_data.domain->blocks[0].patches[0]);

        if(time_interp){
            CHECK(fclaw3dx_clawpatch_get_q_timesync(test_data.glob,&test_data.domain->blocks[0].patches[0],time_interp) == cp->griddata_time_interpolated.dataPtr());
        } else {
            CHECK(fclaw3dx_clawpatch_get_q_timesync(test_data.glob,&test_data.domain->blocks[0].patches[0],time_interp) == cp->griddata.dataPtr());
        }
    }
}
TEST_CASE("fclaw3dx_clawpatch user_data")
{
    SinglePatchDomain test_data;
    test_data.setup();

    //CHECK
    void* user_data = (void *) 1234;
    fclaw3dx_clawpatch_set_user_data(test_data.glob, &test_data.domain->blocks[0].patches[0],user_data);
    CHECK(fclaw3dx_clawpatch_get_user_data(test_data.glob,&test_data.domain->blocks[0].patches[0]) == user_data);
}
TEST_CASE("fclaw3dx_clawpatch solver_data")
{
    SinglePatchDomain test_data;
    test_data.setup();

    //CHECK
    void* user_data = (void *) 1234;
    fclaw3dx_clawpatch_set_solver_data(test_data.glob, &test_data.domain->blocks[0].patches[0],user_data);
    CHECK(fclaw3dx_clawpatch_get_solver_data(test_data.glob,&test_data.domain->blocks[0].patches[0]) == user_data);
}

TEST_CASE("fclaw3dx_clawpatch_size")
{
    for(int mx   : {4,5,6})
    for(int my   : {4,5,6})
    for(int mz   : {4,5,6})
    for(int mbc  : {1,2})
    for(int meqn : {1,2})
    {
        fclaw2d_global_t* glob = fclaw2d_global_new(); 
        fclaw2d_vtables_initialize(glob);

        SinglePatchDomain test_data;
        test_data.opts.mx   = mx;
        test_data.opts.my   = my;
        test_data.opts.mz   = mz;
        test_data.opts.mbc  = mbc;
        test_data.opts.meqn = meqn;
        test_data.setup();

        //CHECK
        CHECK(fclaw3dx_clawpatch_size(test_data.glob) == (size_t) (mx+2*mbc)*(my+2*mbc)*(mz+2*mbc)*meqn);

        fclaw2d_global_destroy(glob);
    }
}

TEST_CASE("fclaw3dx_clawpatch_metric_scalar")
{
    SinglePatchDomain test_data;
    test_data.setup();

    double *mp_area, *mp_edgelengths, *mp_curvature;
    double *area, *edgelengths, *curvature;

    fclaw2d_metric_patch_scalar(test_data.glob, &test_data.domain->blocks[0].patches[0], &mp_area, &mp_edgelengths, &mp_curvature);

    fclaw3dx_clawpatch_metric_scalar(test_data.glob, &test_data.domain->blocks[0].patches[0], &area, &edgelengths, &curvature);

    CHECK(mp_area == area);
    CHECK(mp_edgelengths == edgelengths);
    CHECK(mp_curvature == curvature);
}

TEST_CASE("fclaw3dx_clawpatch_metric_vector")
{
    fclaw2d_global_t* glob = fclaw2d_global_new(); 
    fclaw2d_vtables_initialize(glob);

    SinglePatchDomain test_data;
    test_data.setup();

    double *mp_xnormals, *mp_ynormals, *mp_xtangents, *mp_ytangents, *mp_curvature;
    double *xnormals, *ynormals, *xtangents, *ytangents, *curvature;

    fclaw2d_metric_patch_vector(test_data.glob, &test_data.domain->blocks[0].patches[0], 
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

    fclaw2d_global_destroy(glob);
}

TEST_CASE("fclaw3dx_clawpatch_metric_data")
{
    SinglePatchDomain test_data;
    test_data.setup();

    double *mp_xp, *mp_yp, *mp_zp, *mp_xd, *mp_yd, *mp_zd, *mp_area;
    double *xp, *yp, *zp, *xd, *yd, *zd, *area;

    fclaw2d_metric_patch_mesh_data(test_data.glob, &test_data.domain->blocks[0].patches[0], 
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

TEST_CASE("fclaw3dx_clawpatch_metric_data2")
{
    SinglePatchDomain test_data;
    test_data.setup();

    double *mp_xnormals, *mp_ynormals, *mp_xtangents, *mp_ytangents, *mp_surfnormals, *mp_edgelengths, *mp_curvature;
    double *xnormals, *ynormals, *xtangents, *ytangents, *surfnormals, *edgelengths, *curvature;

    fclaw2d_metric_patch_mesh_data2(test_data.glob, &test_data.domain->blocks[0].patches[0], 
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
namespace{
    double timeinterp_alpha;
    int timeinterp_mint;
    fclaw3dx_clawpatch_t* timeinterp_cp;
}
TEST_CASE("fclaw3dx_clawpatch setup_timeinterp")
{
    SinglePatchDomain test_data;
    test_data.opts.interp_stencil_width=2;
    timeinterp_mint = 2;
    test_data.setup();


    timeinterp_cp = fclaw3dx_clawpatch_get_clawpatch(&test_data.domain->blocks[0].patches[0]);
    timeinterp_alpha = 0.90210;


    fclaw3dx_clawpatch_vtable_t * clawpatch_vt = fclaw3dx_clawpatch_vt(test_data.glob);
    clawpatch_vt->fort_timeinterp = [] (const int *mx, const int *my, const int *mz, 
                                        const int *mbc, const int *meqn, const int *psize, 
                                        double qcurr[], double qlast[], double qinterp[], 
                                        const double *alpha, const int *ierror)
    {
        CHECK(*mx == timeinterp_cp->mx);
        CHECK(*my == timeinterp_cp->my);
        CHECK(*mz == timeinterp_cp->mz);
        CHECK(*mbc == timeinterp_cp->mbc);
        CHECK(*meqn == timeinterp_cp->meqn);
        CHECK(*psize == (*mx)*(*my)*(*mz)
                        - (*mx-2*timeinterp_mint)*(*my-2*timeinterp_mint)*(*mz));
        CHECK(qcurr == timeinterp_cp->griddata.dataPtr());
        CHECK(qlast == timeinterp_cp->griddata_last.dataPtr());
        CHECK(qinterp == timeinterp_cp->griddata_time_interpolated.dataPtr());
        CHECK(*alpha == timeinterp_alpha);
    };

    fclaw2d_patch_setup_timeinterp(test_data.glob, &test_data.domain->blocks[0].patches[0], 
                                   timeinterp_alpha);
}
namespace{
    fclaw3dx_clawpatch_t* t4r_cp;
    int t4r_tag_patch;
    int t4r_init_flag;
}
TEST_CASE("fclaw3dx_clawpatch tag4refinement")
{
    for(int tag_patch : {true,false})
    for(int init_flag : {true,false})
    {
        SinglePatchDomain test_data;
        test_data.fopts.refine_threshold = 0.90210;
        test_data.setup();


        t4r_cp = fclaw3dx_clawpatch_get_clawpatch(&test_data.domain->blocks[0].patches[0]);

        fclaw3dx_clawpatch_vtable_t * clawpatch_vt = fclaw3dx_clawpatch_vt(test_data.glob);

        t4r_tag_patch = tag_patch;
        t4r_init_flag = init_flag;
        clawpatch_vt->fort_tag4refinement = [](const int *mx, const int *my, const int *mz, 
                                               const int *mbc, const int *meqn, 
                                               const double *xlower, const double *ylower, const double *zlower, 
                                               const double *dx, const double *dy, const double *dz, 
                                               const int *blockno, 
                                               double q[], const double *tag_threshold, const int *init_flag, 
                                               int *tag_patch)
        {
            CHECK(*mx == t4r_cp->mx);
            CHECK(*my == t4r_cp->my);
            CHECK(*mz == t4r_cp->mz);
            CHECK(*mbc == t4r_cp->mbc);
            CHECK(*meqn == t4r_cp->meqn);
            CHECK(*xlower == t4r_cp->xlower);
            CHECK(*ylower == t4r_cp->ylower);
            CHECK(*zlower == t4r_cp->zlower);
            CHECK(*dx == t4r_cp->dx);
            CHECK(*dy == t4r_cp->dy);
            CHECK(*dz == t4r_cp->dz);
            CHECK(*blockno == t4r_cp->blockno);
            CHECK(q == t4r_cp->griddata.dataPtr());
            CHECK(*tag_threshold == .90210);
            CHECK(*init_flag == t4r_init_flag);
            *tag_patch = t4r_tag_patch;
        };

        CHECK(fclaw2d_patch_tag4refinement(test_data.glob, &test_data.domain->blocks[0].patches[0], 0, 0, t4r_init_flag) == t4r_tag_patch);
    }
}
TEST_CASE("fclaw3dx_clawpatch tag4refinement negative refine threshold")
{
    for(int init_flag : {true,false})
    {
        SinglePatchDomain test_data;
        test_data.fopts.refine_threshold = -0.90210;
        test_data.setup();

        t4r_init_flag = init_flag;

        CHECK(fclaw2d_patch_tag4refinement(test_data.glob, &test_data.domain->blocks[0].patches[0], 0, 0, t4r_init_flag) == true);
    }
}
namespace{
    fclaw3dx_clawpatch_t* t4c_cp;
    fclaw3dx_clawpatch_t* t4c_cp1;
    fclaw3dx_clawpatch_t* t4c_cp2;
    fclaw3dx_clawpatch_t* t4c_cp3;
    int t4c_tag_patch;
    int t4c_init_flag;
}
TEST_CASE("fclaw3dx_clawpatch tag4coarsening")
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


            t4c_cp = fclaw3dx_clawpatch_get_clawpatch(&test_data.domain->blocks[0].patches[0]);
            t4c_cp1 = fclaw3dx_clawpatch_get_clawpatch(&test_data.domain->blocks[0].patches[1]);
            t4c_cp2 = fclaw3dx_clawpatch_get_clawpatch(&test_data.domain->blocks[0].patches[2]);
            t4c_cp3 = fclaw3dx_clawpatch_get_clawpatch(&test_data.domain->blocks[0].patches[3]);

            fclaw3dx_clawpatch_vtable_t * clawpatch_vt = fclaw3dx_clawpatch_vt(test_data.glob);

            clawpatch_vt->fort_tag4coarsening = [](const int *mx, const int *my, const int *mz, 
                                                   const int *mbc, const int *meqn, 
                                                   double xlower[], double ylower[], double zlower[], 
                                                   const double *dx, const double *dy, const double *dz, 
                                                   const int *blockno, 
                                                   double q0[], double q1[], double q2[], double q3[], 
                                                   const double *tag_threshold, const int *init_flag, int *tag_patch)
            {
                CHECK(*mx == t4c_cp->mx);
                CHECK(*my == t4c_cp->my);
                CHECK(*mz == t4c_cp->mz);
                CHECK(*mbc == t4c_cp->mbc);
                CHECK(*meqn == t4c_cp->meqn);
                CHECK(*xlower == t4c_cp->xlower);
                CHECK(*ylower == t4c_cp->ylower);
                CHECK(*zlower == t4c_cp->zlower);
                CHECK(*dx == t4c_cp->dx);
                CHECK(*dy == t4c_cp->dy);
                CHECK(*dz == t4c_cp->dz);
                CHECK(*blockno == t4c_cp->blockno);
                CHECK(q0 == t4c_cp->griddata.dataPtr());
                CHECK(q1 == t4c_cp1->griddata.dataPtr());
                CHECK(q2 == t4c_cp2->griddata.dataPtr());
                CHECK(q3 == t4c_cp3->griddata.dataPtr());
                CHECK(*tag_threshold == .90210);
                CHECK(*init_flag == t4c_init_flag);
                *tag_patch = t4c_tag_patch;
            };

            CHECK(fclaw2d_patch_tag4coarsening(test_data.glob, &test_data.domain->blocks[0].patches[0], 0, 0, t4c_init_flag) == t4c_tag_patch);
        }
    }
}
TEST_CASE("fclaw3dx_clawpatch tag4coarsening negative coarsen threshold")
{
    for(int init_flag : {true,false})
    {
        t4c_init_flag = init_flag;
    
        QuadDomain test_data;
        test_data.fopts.coarsen_threshold = -0.90210;
        test_data.setup();

        CHECK(fclaw2d_patch_tag4coarsening(test_data.glob, &test_data.domain->blocks[0].patches[0], 0, 0, t4c_init_flag) == false);
    }
}
namespace{
    fclaw3dx_clawpatch_t* i2f_ccp;
    fclaw3dx_clawpatch_t* i2f_cp;
    fclaw3dx_clawpatch_t* i2f_cp1;
    fclaw3dx_clawpatch_t* i2f_cp2;
    fclaw3dx_clawpatch_t* i2f_cp3;
    std::bitset<4> i2f_igrids;
    int i2f_manifold;
}
TEST_CASE("fclaw3dx_clawpatch interpolate2fine")
{
    i2f_manifold = false;

    QuadDomain fine_test_data;
    fine_test_data.fopts.manifold = i2f_manifold;
    fine_test_data.setup();
    SinglePatchDomain coarse_test_data;
    coarse_test_data.fopts.manifold = i2f_manifold;
    coarse_test_data.setup();


    i2f_ccp = fclaw3dx_clawpatch_get_clawpatch(&coarse_test_data.domain->blocks[0].patches[0]);
    i2f_cp = fclaw3dx_clawpatch_get_clawpatch(&fine_test_data.domain->blocks[0].patches[0]);
    i2f_cp1 = fclaw3dx_clawpatch_get_clawpatch(&fine_test_data.domain->blocks[0].patches[1]);
    i2f_cp2 = fclaw3dx_clawpatch_get_clawpatch(&fine_test_data.domain->blocks[0].patches[2]);
    i2f_cp3 = fclaw3dx_clawpatch_get_clawpatch(&fine_test_data.domain->blocks[0].patches[3]);

    fclaw3dx_clawpatch_vtable_t * clawpatch_vt = fclaw3dx_clawpatch_vt(coarse_test_data.glob);

    clawpatch_vt->fort_interpolate2fine = [](const int *mx, const int *my, const int *mz, 
                                             const int *mbc, const int *meqn, 
                                             double qcoarse[], double qfine[], double areacoarse[], double areafine[], 
                                             const int *igrid, const int *manifold)
    {
        CHECK(*mx == i2f_cp->mx);
        CHECK(*my == i2f_cp->my);
        CHECK(*mz == i2f_cp->mz);
        CHECK(*mbc == i2f_cp->mbc);
        CHECK(*meqn == i2f_cp->meqn);

        CHECK(qcoarse == i2f_ccp->griddata.dataPtr());
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

        CHECK(i2f_igrids[*igrid] == false);
        i2f_igrids[*igrid] = true;
        CHECK(*manifold == i2f_manifold);
    };

    fclaw2d_patch_interpolate2fine(coarse_test_data.glob,
                                   &coarse_test_data.domain->blocks[0].patches[0],
                                   &fine_test_data.domain->blocks[0].patches[0],
                                   0, 0, 0);

    CHECK(i2f_igrids.all());
}