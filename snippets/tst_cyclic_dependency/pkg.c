#include "global.h"
#include "pkg.h"

void pkg_add_value(global_t *g, int v)
{
    g->pkg->value = v;
}
