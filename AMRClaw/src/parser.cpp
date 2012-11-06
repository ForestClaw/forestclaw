#include <iostream>
#include <sstream>
#include <vector>

using namespace std;

#include "parser.H"

int     Parser::s_argc = NULL;
char**  Parser::s_argv = NULL;

Parser::Parser()
{
    m_ini = NULL;
}

void Parser::define(int a_argc, char **a_argv)
{
    s_argc = a_argc;
    s_argv = a_argv;
}

Parser::Parser(const char *a_section)
{
    if (s_argv == NULL)
    {
        printf("parser.cpp : parser is not defined\n");
        exit(1);
    }
    // Read second argument, since first argument is the "-F"
    m_ini = iniparser_load(s_argv[2], stderr);
    if (m_ini == NULL)
    {
        fprintf(stderr, "Cannot parse file: %s\n", s_argv[1]);
        exit(1);
    }
    m_section.assign(a_section);
}

Parser::~Parser()
{
    iniparser_freedict(m_ini);
}

void Parser::dump()
{
    iniparser_dump(m_ini, stderr);
}

const char* Parser::append_key(const char* key)
{
    string section_plus_key;
    section_plus_key.assign(m_section);
    section_plus_key.append(":");
    section_plus_key.append(key);
    return section_plus_key.c_str();
}

int Parser::get_int(const char *key, int notfound)
{
    const char* section_plus_key = append_key(key);
    return iniparser_getint(m_ini, section_plus_key, notfound);
}

double Parser::get_double(const char *key, double notfound)
{
    const char* section_plus_key = append_key(key);
    return iniparser_getdouble(m_ini, section_plus_key, notfound);
}

int Parser::get_boolean(const char *key, int notfound)
{
    const char* section_plus_key = append_key(key);
    return iniparser_getboolean(m_ini, section_plus_key, notfound);
}

vector<int> Parser::get_int_array(const char *key)
{

    const char* section_plus_key = append_key(key);
    vector<int> numbers;
    char *line;
    line = iniparser_getstring(m_ini,section_plus_key, NULL);
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
    const char* section_plus_key = append_key(key);
    vector<double> numbers;
    char *line;
    line = iniparser_getstring(m_ini,section_plus_key, NULL);
    if (line == NULL)
    {
        return numbers;
    }

    stringstream lineStream(line);

    double num;
    while (lineStream >> num) numbers.push_back(num);

    return numbers;
}
