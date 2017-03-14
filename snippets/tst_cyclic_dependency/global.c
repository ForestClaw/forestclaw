#include "global.h"
#include "pkg.h"

void global_add_pkg(global_t* g,pkg_t *p)
{
    g->pkg = p;
}

int global_get_value(global_t* g)
{
    return g->pkg->value;
}
