#include <test/catch.hpp>
#include <test/test_config.h>
#include "phasefield_patch_operator.h"
#include <test/utils/DomainReader.h>
#include <ThunderEgg/ValVector.h>
#include <ThunderEgg/BiLinearGhostFiller.h>
#include <ThunderEgg/DomainTools.h>

using namespace std;
using namespace ThunderEgg;

TEST_CASE("default lambda value is 999", "[applications/elliptic/phasefield]")
{
    CHECK(phasefield::getLambda()==999);
}

TEST_CASE("setLambda", "[applications/elliptic/phasefield]")
{
    double orig_lambda = phasefield::getLambda();
    phasefield::setLambda(2);
    CHECK(phasefield::getLambda()==2);
    phasefield::setLambda(orig_lambda);
}
TEST_CASE("constructor sets s array correctly", "[applications/elliptic/phasefield]"){
    double orig_lambda = phasefield::getLambda();
    phasefield::setLambda(-1);

    DomainReader<2> reader(TEST_SRC_DIR "/test/mesh_files/2d_uniform_2x2_mpi1.json", {10,10}, 2);
    auto domain = reader.getCoarserDomain();

    auto phi_n = ValVector<2>::GetNewVector(domain, 3);

    auto ghost_filler = make_shared<BiLinearGhostFiller>(domain, GhostFillingType::Faces);

    fc2d_thunderegg_options_t * mg_opt = new fc2d_thunderegg_options_t();
    mg_opt->boundary_conditions = new int[4];
    mg_opt->boundary_conditions[0] = GENERATE(1,2);
    mg_opt->boundary_conditions[1] = GENERATE(1,2);
    mg_opt->boundary_conditions[2] = GENERATE(1,2);
    mg_opt->boundary_conditions[3] = GENERATE(1,2);


    phasefield_options_t * phase_opt = new phasefield_options_t();

    phasefield op(mg_opt, phase_opt, phi_n, domain, ghost_filler);

    const int* s = op.getS();

    for(int i=0;i<4;i++){
        if(mg_opt->boundary_conditions[i]==1){
            CHECK(s[i]==-1);
        }else{
            CHECK(s[i]==1);
        }
    }

    phasefield::setLambda(orig_lambda);

    //cleanup
    delete[] mg_opt->boundary_conditions;
    delete mg_opt;
    delete phase_opt;
}


double one(std::array<double, 2> coord){
    double x = coord[0];
    double y = coord[1];
    return sin(x)*cos(y);
}
double two(std::array<double, 2> coord){
    double x = coord[0];
    double y = coord[1];
    return cos(x)*sin(y);
}
TEST_CASE("addGhostToRHS", "[applications/elliptic/phasefield]"){
    double orig_lambda = phasefield::getLambda();
    phasefield::setLambda(-1);

    DomainReader<2> reader(TEST_SRC_DIR "/test/mesh_files/2d_uniform_2x2_mpi1.json", {10,10}, 2);
    auto domain = reader.getFinerDomain();

    auto x = ValVector<2>::GetNewVector(domain, 2);
    auto x_just_ghosts = ValVector<2>::GetNewVector(domain, 2);
    auto y_expected = ValVector<2>::GetNewVector(domain, 2);
    auto y = ValVector<2>::GetNewVector(domain, 2);
    auto phi_n = ValVector<2>::GetNewVector(domain, 3);

    auto ghost_filler = make_shared<BiLinearGhostFiller>(domain, GhostFillingType::Faces);

    DomainTools::SetValuesWithGhost<2>(domain,x,one,one);
    DomainTools::SetValues<2>(domain,x_just_ghosts,one,one);
    ghost_filler->fillGhost(x_just_ghosts);
    for(auto pinfo : domain->getPatchInfoVector()){
        auto lds = x_just_ghosts->getLocalDatas(pinfo.local_index);
        for(Side<2> s : Side<2>::getValues()){
            if(pinfo.hasNbr(s)){
                for(int component = 0;component<2;component++){
                    auto inner_slice = lds[component].getSliceOn(s,{0});
                    auto ghost_slice = lds[component].getGhostSliceOn(s,{0});
                    for(int i=0;i<10;i++){
                        ghost_slice[{i}]+=inner_slice[{i}];
                    }
                }
            }
        }
    }
    x_just_ghosts->set(0);

    DomainTools::SetValuesWithGhost<2>(domain,phi_n,two,two,two);

    fc2d_thunderegg_options_t * mg_opt = new fc2d_thunderegg_options_t();
    mg_opt->boundary_conditions = new int[4]{1,1,1,1};

    phasefield_options_t * phase_opt = new phasefield_options_t();
    phase_opt->S = 0.5;
    phase_opt->alpha=400;
    phase_opt->m = 0.035;
    phase_opt->xi = 0.0065;
    phase_opt->k = 4;
    phase_opt->gamma = 0.5;
    phase_opt->r0 = 0.1;

    phasefield op(mg_opt, phase_opt, phi_n, domain, ghost_filler);

    op.apply(x,y);
    for(auto pinfo : domain->getPatchInfoVector()){
        auto x_lds = x_just_ghosts->getLocalDatas(pinfo.local_index);
        auto y_lds = y_expected->getLocalDatas(pinfo.local_index);
        op.applySinglePatch(pinfo,x_lds,y_lds,false);
    }
    y_expected->scaleThenAdd(-1.0,y);

    ghost_filler->fillGhost(x);
    for(auto pinfo : domain->getPatchInfoVector()){
        INFO("x_start: " <<pinfo.starts[0]);
        INFO("y_start: " <<pinfo.starts[1]);
        auto expected_lds = y_expected->getLocalDatas(pinfo.local_index);
        auto lds = y->getLocalDatas(pinfo.local_index);
        auto us = x->getLocalDatas(pinfo.local_index);
        op.addGhostToRHS(pinfo,us,lds);
        for(int component = 0;component<2;component++){
            INFO("Component: "<< component);
            for(int yi = 0; yi<10; yi++){
                INFO("yi: "<< yi);
                for(int xi = 0; xi<10; xi++){
                    INFO("xi: "<< xi);
                    CHECK(lds[component][{xi,yi}]==Approx(expected_lds[component][{xi,yi}]));
                }
            }
        }
    }

    //cleanup
    phasefield::setLambda(orig_lambda);
    delete[] mg_opt->boundary_conditions;
    delete mg_opt;
}