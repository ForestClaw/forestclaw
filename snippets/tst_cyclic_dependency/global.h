#ifndef GLOBAL_H
#define GLOBAL_H

typedef struct global global_t;
struct pkg;

struct global
{
    struct pkg *pkg;
};

void global_add_pkg(global_t* g, struct pkg *p);
int global_get_value(global_t* g);

#endif
