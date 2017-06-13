#include <stdio.h>
#include <stdlib.h>

#include "global.h"
#include "pkg.h"

int main()
{

    global_t s_g, *g = &s_g;
    pkg_t s_pkg, *pkg = &s_pkg;

    pkg->value = 2;

    global_add_pkg(g,pkg);
    int v = global_get_value(g);

    printf("value is %d\n",v);

    pkg_add_value(g,47);
    v = global_get_value(g);

    printf("value is %d\n",v);


    return 0;
}
