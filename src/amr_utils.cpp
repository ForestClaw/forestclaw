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
#include "amr_forestclaw.H"

int pow_int(int a, int n)
{
    int b = a;
    for(int i = 0; i < n-1; i++)
    {
        b *= a;
    }
    return b;
}


global_parms::global_parms()
{
    m_maxmwaves = 20;
    m_mthlim = new int[m_maxmwaves];
    m_refratio = 2;  // hardwired for now.
}

global_parms::~global_parms()
{
    // delete [] m_mthbc;
    delete [] m_mthlim;
}


void global_parms::get_inputParams()
{
    // read data using Fortran file
    int icycle;
    int imanifold;
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
                m_minlevel,
                m_maxlevel,
                icycle,
                imanifold);

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

    // I should just figure out how bool and Fortran 'logical' vars can be passed around.
    m_subcycle = icycle == 1;
    m_manifold = imanifold == 1;
    m_mapped = false;
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
  cout << "subcycling " << (m_subcycle ? "Yes" : "No") << endl;
  cout << "Manifold " << (m_manifold ? "Yes" : "No") << endl;

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

  cout << "Min configured level = " << m_minlevel << endl;
  cout << "Max configured level = " << m_maxlevel << endl;
}


FArrayBox::FArrayBox()
{
    m_data = NULL;
    m_size = 0;
}

FArrayBox::~FArrayBox()
{
    if (m_data != NULL)
    {
        delete [] m_data;
    }
    m_data = NULL;
    m_size = 0;
}

void FArrayBox::define(int a_size,const Box& a_box)
{
    if (a_size == 0)
    {
        printf("FArrayBox::define(int a_size, const Box& a_box \n");
        exit(1);
    }

    if (m_data != NULL)
    {
        delete [] m_data;
    }
    m_data = new Real[a_size];
    m_size = a_size;
    m_box = a_box;
}

void FArrayBox::define(const Box& a_box, int a_fields)
{
    int box_size = (a_box.bigEnd(0) - a_box.smallEnd(0) + 1)*(a_box.bigEnd(1) - a_box.smallEnd(1) + 1);
    define(box_size*a_fields,a_box);
}


// copy constructor
void FArrayBox::operator=(const FArrayBox& fbox)
{
    if (fbox.m_size != m_size)
    {
        if (m_data != NULL)
        {
            delete [] m_data;
            m_data = NULL;
        }
        m_data = new Real[fbox.m_size];
        m_size = fbox.m_size;
    }
    if (m_data == NULL)
    {
        if (m_size == fbox.m_size)
        {
            printf("FArrayBox::operator=() : sizes are equal, but m_data == NULL\n");
            exit(1);
        }
    }
    Real *copy = fbox.m_data;
    m_box = fbox.m_box;
    for (int i = 0; i < fbox.m_size; i++)
    {
        m_data[i] = copy[i];
    }
}


Real* FArrayBox::dataPtr()
{
    return m_data;
}

Box FArrayBox::box()
{
    return m_box;
}

int FArrayBox::size()
{
    return m_size;
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

subcycle_manager::subcycle_manager() {}
subcycle_manager::~subcycle_manager() {}

void subcycle_manager::define(fclaw2d_domain_t *domain,
                              global_parms *gparms,
                              const amr_options_t *amropts,
                              const Real& a_t_curr)
{
    m_t_minlevel = a_t_curr;
    m_refratio = gparms->m_refratio;

    /* query the levels that exist on this processor */
    m_minlevel = domain->minlevel_all;
    m_maxlevel = domain->maxlevel_all;
    m_levels.resize(m_maxlevel + 1);

    bool subcycle = gparms->m_subcycle;
    m_nosubcycle = !subcycle;
    for (int level = m_minlevel; level <= m_maxlevel; level++)
    {
        int np = num_patches(domain,level);
        m_levels[level].define(level,m_refratio,np,m_maxlevel,a_t_curr,subcycle);
    }
}

bool subcycle_manager::nosubcycle()
{
    return m_nosubcycle;
}

void subcycle_manager::set_dt_minlevel(const Real& a_dt_minlevel)
{
    // Time step for minimum level (i.e. coarsest non-empty level).
    m_dt_minlevel = a_dt_minlevel;

    Real dt_level = a_dt_minlevel;
    m_levels[m_minlevel].set_dt(dt_level);
    for (int level = m_minlevel+1; level <= m_maxlevel; level++)
    {
        if (!nosubcycle())
        {
            dt_level /= m_refratio;
        }
        m_levels[level].set_dt(dt_level);
    }
}


int subcycle_manager::minlevel_factor()
{
    int factor = 1;
    for (int level = 1; level <= m_minlevel; level++)
    {
        factor *= m_refratio;
    }
    return factor;
}

int subcycle_manager::maxlevel_factor()
{
    int factor = 1;
    for (int level = 1; level <= m_maxlevel; level++)
    {
        factor *= m_refratio;
    }
    return factor;
}


int subcycle_manager::minlevel()
{
    return m_minlevel;
}

int subcycle_manager::maxlevel()
{
    return m_maxlevel;
}

int subcycle_manager::last_step(const int& a_level)
{
    return m_levels[a_level].m_last_step;
}

void subcycle_manager::increment_step_counter(const int& a_level)
{
    m_levels[a_level].increment_step_counter();
}

void subcycle_manager::increment_time(const int& a_level)
{
    m_levels[a_level].increment_time();
}

bool subcycle_manager::is_coarsest(const int& a_level)
{
    return a_level == m_minlevel;
}

bool subcycle_manager::is_finest(const int& a_level)
{
    return a_level == m_maxlevel;
}

Real subcycle_manager::dt(const int& a_level)
{
    return m_levels[a_level].dt();
}

bool subcycle_manager::can_advance(const int& a_level, const int& a_curr_step)
{
    bool verbose = false;
    bool b1 = solution_updated(a_level,a_curr_step); // do we need this?  We shouldn't be here if
                                                     // we have not taken a time step to 'a_curr_step'
    bool b2 = level_exchange_done(a_level);
    bool b3 = exchanged_with_coarser(a_level);
    bool b4 = exchanged_with_finer(a_level);  // This may not be needed.
    if (verbose)
    {
        if (!b1)
        {
            cout << " --> (can_advance) Solution at level " << a_level <<
                " has not been updated at step " << a_curr_step << endl;
        }
        if (!b2)
        {
            cout << " --> (can_advance) Level exchange at level " << a_level <<
                " has not been done at step " << a_curr_step << endl;
        }
        if (!b3)
        {
            cout << " --> (can_advance) Exchange with coarse grid at level " << a_level-1 <<
                " not done at step " << a_curr_step << endl;
        }
        if (!b4)
        {
            cout << " --> (can_advance) Exchange with finer grid at level " << a_level+1 <<
                " not done at step " << a_curr_step << endl;
        }
    }
    return b1 && b2 && b3 && b4;
}

Real subcycle_manager::current_time(const int& a_level)
{
    return m_levels[a_level].current_time();
}

bool subcycle_manager::solution_updated(const int& a_level, const int& a_step)
{
    return m_levels[a_level].m_last_step >= a_step;
}

bool subcycle_manager::level_exchange_done(const int& a_level)
{
    return m_levels[a_level].level_exchange_done();
}

bool subcycle_manager::exchanged_with_coarser(const int& a_level)
{
    if (is_coarsest(a_level))
        return true;
    else
        return m_levels[a_level].exchanged_with_coarser();
}

bool subcycle_manager::exchanged_with_finer(const int& a_level)
{
    if (is_finest(a_level))
        return true;
    else
        return m_levels[a_level].exchanged_with_finer();
}

void subcycle_manager::increment_level_exchange_counter(const int& a_level)
{
    m_levels[a_level].increment_level_exchange_counter();
}

void subcycle_manager::increment_coarse_exchange_counter(const int& a_level)
{
    if (!is_coarsest(a_level))
        m_levels[a_level].increment_coarse_exchange_counter();
}

void subcycle_manager::increment_fine_exchange_counter(const int& a_level)
{
    if (!is_finest(a_level))
        m_levels[a_level].increment_fine_exchange_counter();
}

int subcycle_manager::step_inc(const int& a_level)
{
    return m_levels[a_level].m_step_inc;
}



level_data::level_data() { }
level_data::~level_data() {}

void level_data::define(const int& a_level,
                        const int& a_refratio,
                        const int& a_patches_at_level,
                        const int& a_maxlevel,
                        const Real& a_time,
                        const bool& a_subcycle)
{
    m_level = a_level;
    m_last_step = 0;
    m_last_level_exchange = 0;   // Assume that the level exchange has been done at subcycled time
                                 // step 0.
    m_num_patches = a_patches_at_level;
    m_time = a_time;


    // This factor determinines how many finest level steps are equivalent to a single step at
    // this level.
    // Example : 2 level 1 steps are equal to 1 level 0 step, so we have
    //         m_step_inc = 2 for level 0
    //         m_step_inc = 1 for level 1
    m_step_inc = 1;
    if (a_subcycle)
    {
        for (int level = a_maxlevel; level > a_level; level--)
        {
            m_step_inc *= a_refratio;
        }
    }
    m_last_coarse_exchange = -m_step_inc;
    m_last_fine_exchange = -m_step_inc;

}

void level_data::set_dt(const Real& a_dt_level)
{
    m_dt = a_dt_level;
}

Real level_data::dt()
{
    return m_dt;
}

void level_data::increment_step_counter()
{
    m_last_step += m_step_inc;
}

void level_data::increment_time()
{
    m_time += m_dt;
}

Real level_data::current_time()
{
    return m_time;
}

void level_data::increment_level_exchange_counter()
{
    m_last_level_exchange += m_step_inc;
}

void level_data::increment_coarse_exchange_counter()
{
    m_last_coarse_exchange += m_step_inc;
}

void level_data::increment_fine_exchange_counter()
{
    m_last_fine_exchange += m_step_inc;
}


bool level_data::level_exchange_done()
{
    return m_last_level_exchange == m_last_step;
}

bool level_data::exchanged_with_coarser()
{
    return m_last_coarse_exchange == m_last_step;
}

bool level_data::exchanged_with_finer()
{
    return m_last_fine_exchange == m_last_step;
}


void cb_num_patches(fclaw2d_domain_t *domain,
	fclaw2d_patch_t *patch, int block_no, int patch_no, void *user)
{
  (*(int *) user)++;
}

int num_patches(fclaw2d_domain_t *domain, int level)
{
  int count = 0;
  fclaw2d_domain_iterate_level(domain, level,
                               cb_num_patches,
                               &count);
  return count;
}
