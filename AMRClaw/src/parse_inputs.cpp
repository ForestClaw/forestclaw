#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#include <iniparser.h>

#include "amr_options.h"
#include "parser.H"

using namespace std;

/* Move this to swirl.cpp, etc */
/*
int main(int argc, char * argv[])
{
    int status;
    // Set static variables
    Parser P;
    P.define(argc,argv);
    status = parse_ini_file();
    return status;
}
*/

/* Read in vector valued inputs */
int parse_ini_file(amr_options_t *amropt)
{
    // Open parser for file stored in s_argc.
    Parser P("Options");

    vector<int> order;
    order = P.get_int_array("order");
    if (order.size() != 2)
    {
        printf("Wrong number of values for vector \"order\"\n");
        exit(1);
    }
    for (int j = 0; j < 2; j++)
    {
        amropt->order[j] = order[j];
    }

    vector<int> mthbc;
    mthbc = P.get_int_array("mthbc");
    if ((int) mthbc.size() != CubeFaces)
    {
        printf("Wrong number of values for vector \"mthbc\"\n");
        exit(1);
    }
    for (int j = 0; j < (int) mthbc.size(); j++)
    {
        amropt->mthbc[j] = mthbc[j];
    }

    vector<int> mthlim;
    mthlim = P.get_int_array("mthlim");
    if ((int) mthlim.size() != amropt->mwaves)
    {
        printf("Wrong number of limiters read in; mwaves = %d\n",amropt->mwaves);
        exit(1);
    }
    amropt->mthlim = new int[amropt->mwaves];
    for (int j = 0; j < amropt->mwaves; j++)
    {
        amropt->mthlim[j] = mthlim[j];
    }

    // Set up 'method' vector used by Clawpack.
    amropt->method[0] = 0; // not used (yet) in forestclaw

    amropt->method[1] = amropt->order[0];
    if (SpaceDim == 2)
    {
        amropt->method[2] = amropt->order[1];
    }
    else
    {
        amropt->method[2] = 10*amropt->order[1] + amropt->order[2];
    }
    amropt->method[3] = amropt->verbosity;
    amropt->method[4] = amropt->src_term;
    amropt->method[5] = amropt->mcapa;
    amropt->method[6] = amropt->maux;

    return 0;
}
