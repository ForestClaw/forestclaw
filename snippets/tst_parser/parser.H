#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#include <iniparser.h>

#include <iostream>
#include <sstream>
#include <vector>

using namespace std;

class Parser
{
public:
    Parser(const char *ini_name);
    ~Parser();
    void dump();

    int get_int(const char *key, int notfound);
    double get_double(const char *key, double notfound);
    int get_boolean(const char *key, int notfound);

    vector<int> get_int_array(const char *key);
    vector<double> get_double_array(const char *key);

protected:
    dictionary	*m_ini;
};
