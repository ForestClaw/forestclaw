#include <test/catch.hpp>
#include <applications/elliptic/phasefield/phasefield_patch_operator.h>
#include <test/utils/DomainReader.h>
#include <ThunderEgg/ValVector.h>
#include <ThunderEgg/BiLinearGhostFiller.h>

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

    DomainReader<2> reader(DATADIR "/test/mesh_files/2d_uniform_2x2_mpi1.json", {10,10}, 2);
    auto domain = reader.getCoarserDomain();

    auto phi_n = ValVector<2>::GetNewVector(domain, 3);

    auto ghost_filler = make_shared<BiLinearGhostFiller>(domain);

    fclaw2d_global_t * glob = fclaw2d_global_new();

    fc2d_thunderegg_options_t * mg_opt = new fc2d_thunderegg_options_t();
    mg_opt->boundary_conditions = new int[4];
    mg_opt->boundary_conditions[0] = GENERATE(1,2);
    mg_opt->boundary_conditions[1] = GENERATE(1,2);
    mg_opt->boundary_conditions[2] = GENERATE(1,2);
    mg_opt->boundary_conditions[3] = GENERATE(1,2);


    phasefield_options_t * phasefield_options = new phasefield_options_t();

    phasefield op(mg_opt, phasefield_options, phi_n, domain, ghost_filler);

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
    fclaw2d_global_destroy(glob);
}