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

global_parms::global_parms()
{
    m_maxmwaves = 20;
    m_mthlim = new int[m_maxmwaves];
    m_mthbc = new int[4];
    m_order = new int[2];
    m_method = new int[7];
}

global_parms::~global_parms()
{
    delete [] m_mthlim;
    delete [] m_mthbc;
    delete [] m_order;
    delete [] m_method;
}


/*
void global_parms::set_order(int *a_order)
{
    m_method[0] = 0; // not used in forestclaw
    m_method[1] = a_order[0];
    if (SpaceDim == 2)
    {
        m_method[2] = a_order[1];
    }
    else
    {
        m_method[2] = 10*a_order[1] + a_order[2];
    }
}

void global_parms::set_src(int a_src)
{
  m_method[4] = a_src;
}

void global_parms::set_mcapa(int a_mcapa)
{
    m_mcapa = a_mcapa;
    m_method[5] = a_mcapa;
}

void global_parms::set_maux(int a_maux)
{
    m_maux = a_maux;
    m_method[6] = a_maux;
}
*/


void global_parms::get_inputParams()
{

    // read data using Fortran file
    inputparms_(m_mx_leaf,
                m_my_leaf,
                m_initial_dt,
                m_tfinal,
                m_max_cfl,
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
                m_order);

    // set_order(m_order);
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

    // set_src(m_src_term);
    // set_mcapa(m_mcapa);
    // set_maux(m_maux);

}


void global_parms::print_inputParams()
{
  cout << endl;
  cout << "CLAWPACK PARAMETERS : " << endl;
  cout << "Initial dt " << m_initial_dt << endl;
  cout << "maximum cfl " << m_max_cfl << endl;
  // cout << "method(1) (fixed time step)       = " << fixedDt << endl;
  cout << "method(2:3) (order of integration) = " << m_method[1] << " " << m_method[2] << endl;
  cout << "method(4) (verbosity) = " << m_method[3] << endl;
  cout << "method(5) (source term splitting) = " << m_method[4] << endl;
  cout << "method(6) (mcapa) = " << m_method[5] << endl;
  cout << "method(7) (maux) = " << m_method[6] << endl;
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

  // cout << "Domain values " << endl;
  // cout << "xlower, xupper = " << m_xlower << " " << m_xupper << endl;
  // cout << "ylower, yupper = " << m_ylower << " " << m_yupper << endl;

#if CH_SPACEDIM == 3
  // cout << "zlower, zupper = " << m_zlower << " " << m_zupper << endl;
#endif

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


}
