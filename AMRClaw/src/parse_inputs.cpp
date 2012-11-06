#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#include <iniparser.h>

#include "amr_options.h"
#include "parser.H"

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

#ifdef __cplusplus
extern "C"
{
#if 0
}                               /* need this because indent is dumb */
#endif
#endif


/* Read in vector valued inputs */
int parse_ini_file(amr_options_t *amropt)
{
    // Open parser for file stored in s_argc.
    Parser P("Options");

    vector<int> order;
    order = P.get_int_array("order");
    if (order.size() != 2)
    {
        printf("Wrong number of values for vector ""order""\n");
        exit(1);
    }
    for (int j = 0; j < 2; j++)
    {
        amropt->order[j] = order[j];
    }

    vector<int> mthbc;
    order = P.get_int_array("mthbc");
    if (order.size() != 4)
    {
        printf("Wrong number of values for vector ""mthbc""\n");
        exit(1);
    }
    for (int j = 0; j < mthbc.size(); j++)
    {
        amropt->mthbc[j] = mthbc[j];
    }

    vector<int> mthlim;
    mthlim = P.get_int_array("mthlim");
    if (mthlim.size() != amropt->mwaves)
    {
        printf("Wrong number of limiters read in; mwaves = %d\n",amropt->mwaves);
        exit(1);
    }
    amropt->mthlim = new int[amropt->mwaves];
    for (int j = 0; j < amropt->mwaves; j++)
    {
        amropt->mthlim[j] = mthlim[j];
    }

    return 0;
}

#ifdef __cplusplus
#if 0
{                               /* need this because indent is dumb */
#endif
}
#endif
