/*
Copyright (c) 2012 Carsten Burstedde, Donna Calhoun
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

 * Redistributions of source code must retain the above copyright notice, this
list of conditions and the following disclaimer.
 * Redistributions in binary form must reproduce the above copyright notice,
this list of conditions and the following disclaimer in the documentation
and/or other materials provided with the distribution.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED.  IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

#ifndef AMR_UTILS_H
#define AMR_UTILS_H

#include <iostream>
#include <cstdlib>
#include <vector>

using std::ifstream;
using std::ios;
using std::vector;

#define CH_assert(a)

#define CH_SPACEDIM 2

#define SpaceDim 2

// Needs to be fixed!
#define Max(a,b)  a

#define cout std::cout

#define endl std::endl

typedef double Real;


extern "C"
{
    void inputparms_(int& mx_leaf,
                     int& my_leaf,
                     Real& initial_dt,
                     Real& tfinal,
                     Real& max_cfl,
                     int& nout,
                     int& src_term,
                     int& verbose,
                     int& mcapa,
                     int& maux,
                     int& meqn,
                     int& mwaves,
                     int& maxmwaves,
                     int  mthlim[],
                     int& mbc,
                     int  mthbc[],
                     int  order[]);
}

class global_parms
{
public:
    int m_mx_leaf;
    int m_my_leaf;

    double m_initial_dt;
    double m_tfinal;
    double m_max_cfl;
    int m_nout;
    int m_src_term;
    int m_verbose;
    int m_mcapa;
    int m_maux;
    int m_meqn;
    int m_mwaves;
    int m_maxmwaves;
    int *m_mthlim;
    int m_mbc;
    int *m_mthbc;
    int *m_order;
    int *m_method;

    global_parms();
    ~global_parms();
    // void set_order(int *a_order);
    // void set_src(int a_src);
    // void set_mcapa(int a_mcapa);
    // void set_maux(int a_maux);

    void get_inputParams(); // Read from input file
    void print_inputParams();  // Print out for beginning info

};

class Box
{
public:
    int smallEnd(int idir) const
        {
            return m_ll[idir];
        }

    int bigEnd(int idir) const
        {
            return m_ur[idir];
        }

private:
    int m_ll[2];
    int m_ur[2];
};

class FArrayBox
{
public:
    Real* dataPtr() const
        {
            return m_data;
        }

    Box box()
        {
            return m_box;
        }

private:
    Real *m_data;
    Box m_box;
};

class ProblemDomain
{
public:
    bool isPeriodic(int idir) const
        {
            return false;
        }

};

#endif
