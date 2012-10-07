#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#include <iniparser.h>

static int  parse_ini_file();

#include "parser.H"

int main(int argc, char * argv[])
{
    int status;
    // Set static variables
    Parser P;
    P.define(argc,argv);
    status = parse_ini_file();
    return status;
}

int parse_ini_file()
{
    // Open parser for section "fclaw"
    Parser P("fclaw");

    P.dump();

    // Get fclaw values
    printf("\n");
    printf("[fclaw]\n");

    int mx;
    mx = P.get_int("mx",0);
    printf("%15s    [%d]\n", "mx", mx);

    int my;
    my = P.get_int("my", 0);
    printf("%15s    [%d]\n", "my", my);

    double tfinal;
    tfinal = P.get_double("tfinal", 0.0);
    printf("%15s    [%g]\n", "Tfinal", tfinal);

    int sub;
    sub = P.get_boolean("subcycling",-1);
    printf("%15s    [%d]\n","Subcycling", sub);

    int mwaves;
    mwaves = P.get_int("mwaves",0);
    printf("%15s    [%d]\n","mwaves",mwaves);

    printf("\n");
    vector<int> mthlim;
    mthlim = P.get_int_array("mthlim");
    if (mthlim.size() != mwaves)
    {
        printf("Wrong number of limiters read in; mwaves = %d\n",mwaves);
        exit(1);
    }
    for (int j = 0; j < mthlim.size(); j++)
    {
        printf("mthlim[%d]   [%d]\n",j,mthlim[j]);
    }

    printf("\n");
    vector<int> mthbc;
    mthbc = P.get_int_array("mthbc");
    if (mthbc.size() != 4)
    {
        printf("Wrong number of boundary conditions\n");
        exit(1);
    }
    for (int j = 0; j < mthbc.size(); j++)
    {
        printf("mthbc[%d]    [%d]\n",j,mthbc[j]);
    }
    printf("etc....\n\n");


    Parser P_user("User");
    printf("[User]\n");
    vector<double> darray;
    darray = P_user.get_double_array("darray");
    for (int j = 0; j < darray.size(); j++)
    {
        printf("darray[%d]    [%g]\n",j,darray[j]);
    }

    printf("\n");
    printf("and what happens when items are missing ...\n");
    int N;
    int badval = 47;
    N = P_user.get_int("N",badval);
    if (N == badval)
    {
        printf("N not found\n");
    }

    vector<int> v;
    v = P_user.get_int_array("v");
    if (v.size() == 0)
    {
        printf("Vector v not found or has length 0.\n");
    }

    vector<int> w;
    w = P_user.get_int_array("w");
    if (w.size() == 0)
    {
        printf("Vector w not found or has length 0.\n");
    }

    return 0;
}
