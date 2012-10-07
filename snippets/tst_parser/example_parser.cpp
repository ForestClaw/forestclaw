#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#include <iniparser.h>

static int  parse_ini_file(const char * ini_name);

#include "parser.H"

int main(int argc, char * argv[])
{
    int status;

    status = parse_ini_file(argv[1]);
    return status;
}

int parse_ini_file(const char * ini_name)
{
    Parser P(ini_name);

    P.dump();

    // Get fclaw values
    printf("\n");
    printf("[fclaw]\n");

    int mx;
    mx = P.get_int("fclaw:mx",0);
    printf("%-15s [%d]\n", "mx:", mx);

    int my;
    my = P.get_int("fclaw:my", 0);
    printf("%-15s [%d]\n", "my:", my);

    double tfinal;
    tfinal = P.get_double("fclaw:tfinal", 0.0);
    printf("%-15s [%g]\n", "Tfinal:", tfinal);

    int sub;
    sub = P.get_boolean("fclaw:subcycling",-1);
    printf("%-15s [%d]\n","Subcycling:", sub);

    int mwaves;
    mwaves = P.get_int("fclaw:mwaves",0);
    printf("%-15s [%d]\n","mwaves",mwaves);

    vector<int> mthlim;
    mthlim = P.get_int_array("fclaw:mthlim");
    if (mthlim.size() != mwaves)
    {
        printf("Wrong number of limiters read in; mwaves = %d\n",mwaves);
        exit(1);
    }
    for (int j = 0; j < mthlim.size(); j++)
    {
        printf("mthlim[%d]:      [%d]\n",j,mthlim[j]);
    }
    printf("etc....\n");

    printf("\n");
    printf("[User]\n");
    vector<double> darray;
    darray = P.get_double_array("User:darray");
    for (int j = 0; j < darray.size(); j++)
    {
        printf("darray[%d]: [%g]\n",j,darray[j]);
    }

    printf("\n");
    printf("and what happens when items are missing ...\n");
    int N;
    int badval = 47;
    N = P.get_int("User:N",badval);
    if (N == badval)
    {
        printf("N not found\n");
    }

    vector<int> v;
    v = P.get_int_array("User:v");
    if (v.size() == 0)
    {
        printf("Vector v not found or has length 0.\n");
    }

    vector<int> w;
    w = P.get_int_array("User:w");
    if (w.size() == 0)
    {
        printf("Vector w not found or has length 0.\n");
    }

    return 0;
}
