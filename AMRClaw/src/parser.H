#ifndef PARSER_H
#define PARSER_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#include <iniparser.h>

#include <iostream>
#include <sstream>
#include <vector>

#include "amr_options.h"

#ifdef __cplusplus
extern "C"
{
#if 0
}                               /* need this because indent is dumb */
#endif
#endif

int parse_ini_options(amr_options_t *amropt);
int parse_ini_output(amr_options_t *amropt,int outstyle);

#ifdef __cplusplus
#if 0
{                               /* need this because indent is dumb */
#endif
}
#endif

class Parser
{
public:
    Parser();
    Parser(const char *section);
    ~Parser();
    void define(int a_argc, char **a_argv);
    void dump();

    int get_int(const char *key, int notfound);
    double get_double(const char *key, double notfound);
    int get_boolean(const char *key, int notfound);

    std::vector<int> get_int_array(const char *key);
    std::vector<double> get_double_array(const char *key);

protected:
    dictionary	*m_ini;
    static int    s_argc;
    static char** s_argv;
    std::string m_section;
private:
    const char* append_key(const char* key);
};

#endif
