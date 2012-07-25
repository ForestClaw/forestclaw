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

#include "amr_utils.H"
#include "forestclaw2d.h"


global_parms::global_parms()
{
    m_maxmwaves = 20;
    m_mthlim = new int[m_maxmwaves];
    m_mthbc = new int[2*SpaceDim];
    m_refratio = 2;  // hardwired for now.
}

global_parms::~global_parms()
{
    delete [] m_mthbc;
    delete [] m_mthlim;
}


void global_parms::get_inputParams()
{
    // read data using Fortran file
    inputparms_(m_mx_leaf,
                m_my_leaf,
                m_initial_dt,
                m_tfinal,
                m_max_cfl,
                m_desired_cfl,
                m_nout,
                m_src_term,
                m_verbose,
                m_mcapa,
                m_maux,
                m_meqn,
                m_mwaves,
                m_maxmwaves,
                m_mthlim,
                m_mbc,
                m_mthbc,
                m_order,
                m_maxlevel);

    // Check maxlevel :
    if (m_maxlevel > P4EST_MAXLEVEL)
    {
        cout << "get_inputParms (amr_utils.f) : User 'maxlevel' > P4EST_MAXLEVEL" << endl;
        exit(1);
    }

    cout << "here... 1" << endl;

    // Set up arrays needed by clawpack.
    m_method[0] = 0; // not used in forestclaw

    m_method[1] = m_order[0];
    if (SpaceDim == 2)
    {
        m_method[2] = m_order[1];
    }
    else
    {
        m_method[2] = 10*m_order[1] + m_order[2];
    }
    m_method[3] = m_verbose;
    m_method[4] = m_src_term;
    m_method[5] = m_mcapa;
    m_method[6] = m_maux;
}


void global_parms::print_inputParams()
{
  cout << endl;
  cout << "CLAWPACK PARAMETERS : " << endl;
  cout << "Initial dt " << m_initial_dt << endl;
  cout << "maximum cfl " << m_max_cfl << endl;
  cout << "desired cfl " << m_desired_cfl << endl;
  cout << "method(2:3) (order of integration) = " << m_method[1] << " " << m_method[2] << endl;
  cout << "method(4) (verbosity) = " << m_method[3] << endl;
  cout << "method(5) (source term splitting) = " << m_method[4] << endl;
  cout << "method(6) (mcapa) = " << m_method[5] << endl;
  cout << "method(7) (maux) = " << m_method[6] << endl;
  cout << endl;

  cout << "refratio (fixed) = " << m_refratio << endl;
  cout << endl;

  cout << "meqn (number of equations) = " << m_meqn << endl;
  cout << "maux (number of auxiliary variables) = " << m_maux << endl;
  cout << "mcapa (location of capacity function) = " << m_mcapa << endl;
  cout << "mwaves (number of waves) = " << m_mwaves << endl;

  cout << "mthlim(mwaves) (limiters) = ";
  for(int i = 0; i < m_mwaves; i++)
    {
      cout << m_mthlim[i] << " ";
    }
  cout << endl << endl;

  cout << "mbc (number of ghost cells) = " << m_mbc << endl;

  /*
  cout << "Auxiliary array type : " << endl;
  for (int i = 0; i < m_maux; i++)
    {
      cout << "  " << m_auxtype[i] << endl;
    }
  cout << endl;
  */


  cout << "mthbc(2*dim) (boundary conditions) = ";
  for(int i = 0; i < 2*SpaceDim; i++)
    {
      cout << m_mthbc[i] << " ";
    }
  cout << endl << endl;

  cout << "Max level = " << m_maxlevel << endl;

}


FArrayBox::FArrayBox()
{
    m_data = NULL;
}

FArrayBox::~FArrayBox()
{
    if (m_data != NULL) {delete [] m_data; m_data = NULL;}
}

void FArrayBox::define(int a_size,const Box& a_box)
{
    m_data = new double[a_size];
    m_box = a_box;
}

Real* FArrayBox::dataPtr() const
{
    return m_data;
}

Box FArrayBox::box()
{
    return m_box;
}

Box::Box()
{
}

Box::Box(const Box& a_box)
{
    for(int idir = 0; idir < 2; idir++)
    {
        m_ll[idir] = a_box.m_ll[idir];
        m_ur[idir] = a_box.m_ur[idir];
    }
}

Box::Box(const int ll[], const int ur[])
{
    for(int idir = 0; idir < 2; idir++)
    {
        m_ll[idir] = ll[idir];
        m_ur[idir] = ur[idir];
    }
}

int Box::smallEnd(int idir) const
{
    return m_ll[idir];
}

int Box::bigEnd(int idir) const
{
    return m_ur[idir];
}
