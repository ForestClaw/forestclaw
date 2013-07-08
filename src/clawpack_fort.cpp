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


/* Difference in nan values :
   The first one is not trapped; the second one is.

   (gdb) print qnan
   $13 = nan(0x8000000000000)
   (gdb) print snan
   $14 = nan(0x000000001)
*/

static
void set_qnan(double& f)
{
    /*
     The quiet nan from math.h
    */
    f = NAN;
}


static
void set_snan(double& f)
{
    /* From :
      "NaNs, Uninitialized Variables, and C++"
      http://codingcastles.blogspot.fr/2008/12/nans-in-c.html
    */
    *((long long*)&f) = 0x7ff0000000000001LL;
}

FArrayBox::FArrayBox()
{
    m_data = NULL;
    m_box = Box();
    m_fields = 0;
    m_size = 0;
}

FArrayBox::~FArrayBox()
{
    if (m_data != NULL)
    {
        delete [] m_data;
        m_data = NULL;
    }
}

void FArrayBox::set_dataPtr(int a_size)
{
    if (a_size < 0)
    {
        printf("FArrayBox::set_dataPtr() : size < 0\n");
        exit(1);
    }
    else if (a_size == 0)
    {
        if (m_data != NULL)
        {
            delete [] m_data;
        }
        m_data = NULL;
    }
    else
    {
        if (m_size != a_size)
        {
            delete [] m_data;
            m_data = new double[a_size];
        }
        else
        {
            // Don't do anything;  m_data is already the right size
        }
        double qnan, snan;
        set_qnan(qnan);
        set_snan(snan);
        for(int i = 0; i < a_size; i++)
        {
            m_data[i] = snan;
        }
    }
}

void FArrayBox::define(const Box& a_box, int a_fields)
{
    int box_size = (a_box.bigEnd(0) - a_box.smallEnd(0) + 1)*
                   (a_box.bigEnd(1) - a_box.smallEnd(1) + 1);
    if (box_size < 0)
    {
        printf("FArrayBox::define(a_box,a_fields) : size < 0\n");
        printf("Box::ll = (%d %d)\n",a_box.smallEnd(0), a_box.smallEnd(1));
        printf("Box::ur = (%d %d)\n",a_box.bigEnd(0), a_box.bigEnd(1));
        exit(1);
    }

    set_dataPtr(box_size*a_fields);

    m_size = box_size*a_fields;
    m_fields = a_fields;
    m_box = a_box;
}

void FArrayBox::copyToMemory(double *data)
{
    for(int i = 0; i < m_size; i++)
    {
        data[i] = m_data[i];
    }
}

void FArrayBox::copyFromMemory(double *data)
{
    for(int i = 0; i < m_size; i++)
    {
        m_data[i] = data[i];
    }
}


// copy constructor
void FArrayBox::operator=(const FArrayBox& fbox)
{
    // reset pointer, if necessary
    set_dataPtr(fbox.m_size);

    m_box = fbox.m_box;
    m_size = fbox.m_size;
    m_fields = fbox.m_fields;

    // copy data to this pointer.
    for (int i = 0; i < fbox.m_size; i++)
    {
        m_data[i] = fbox.m_data[i];
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

int FArrayBox::fields()
{
    return m_fields;
}

// Rename to "IndexBox"
Box::Box()
{
    for (int idir = 0; idir < SpaceDim; idir++)
    {
        m_ll[idir] = 0;
        m_ur[idir] = 0;
    }
}

Box::Box(const Box& a_box)
{
    for(int idir = 0; idir < SpaceDim; idir++)
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
