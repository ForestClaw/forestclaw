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

#include "clawpack_fort.H"

void set_problem_parameters()
{
    setprob_();
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
    m_data = new double[a_size];
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
        m_data = new double[fbox.m_size];
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
    double *copy = fbox.m_data;
    m_box = fbox.m_box;
    for (int i = 0; i < fbox.m_size; i++)
    {
        m_data[i] = copy[i];
    }
}


double* FArrayBox::dataPtr()
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
