#include <iostream>
#include <sstream>
#include <vector>

using namespace std;

#include "parser.H"

Parser::Parser(const char *ini_name)
{
    m_ini = iniparser_load(ini_name, stderr);
    if (m_ini == NULL)
    {
        fprintf(stderr, "Cannot parse file: %s\n", ini_name);
        exit(1);
    }
}

Parser::~Parser()
{
    iniparser_freedict(m_ini);
}

void Parser::dump()
{
    iniparser_dump(m_ini, stderr);
}


int Parser::get_int(const char *key, int notfound)
{
    return iniparser_getint(m_ini, key, notfound);
}

double Parser::get_double(const char *key, double notfound)
{
    return iniparser_getdouble(m_ini, key, notfound);
}

int Parser::get_boolean(const char *key, int notfound)
{
    return iniparser_getboolean(m_ini, key, notfound);
}

vector<int> Parser::get_int_array(const char *key)
{
    vector<int> numbers;
    char *line;
    line = iniparser_getstring(m_ini,key, NULL);
    if (line == NULL)
    {
        return numbers;
    }

    stringstream lineStream(line);

    int num;
    while (lineStream >> num) numbers.push_back(num);

    return numbers;
}

vector<double> Parser::get_double_array(const char *key)
{
    vector<double> numbers;
    char *line;
    line = iniparser_getstring(m_ini,key, NULL);
    if (line == NULL)
    {
        return numbers;
    }

    stringstream lineStream(line);

    double num;
    while (lineStream >> num) numbers.push_back(num);

    return numbers;
}
