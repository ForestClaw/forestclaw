#ifndef PKG_H
#define PKG_H

typedef struct pkg pkg_t;
struct global;

struct pkg
{
    int value;
};

void pkg_add_value(struct global *g, int v);

#endif
