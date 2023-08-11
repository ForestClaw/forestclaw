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

#include <fclaw2d_defs.h>
#include <sc.h>

#include <fclaw_farraybox.hpp>

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

static
void set_big_number(double& f)
{
    f = 1e99;
}


void fclaw_farraybox_set_to_nan(double& f)
{
    set_snan(f);
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

FArrayBox::FArrayBox(const FArrayBox& A)
{
    m_box = A.m_box;
    m_size = A.m_size;
    m_fields = A.m_fields;
    /* A re-allocation, but should hopefully, memory management should be okay */
    set_dataPtr(A.m_size);
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
            if (m_data != NULL)
            {
                delete [] m_data;
                m_data = NULL;
            }
            m_data = new double[a_size];
        }
        else
        {
            // Don't do anything;  m_data is already the right size
        }
        double qnan, snan, big_number;
        // qnan = quiet nan
        // snan = signaling nan
        set_qnan(qnan);
        set_snan(snan);
        set_big_number(big_number);
#if FCLAW_DEBUG
        for(int i = 0; i < a_size; i++)
        {
            m_data[i] = snan;
        }
#endif
    }
}

void FArrayBox::define(const Box& a_box, int a_fields)
{
    int box_size = 1;

    for (int i = 0; i < a_box.boxDim(); ++i)
    {
        box_size *= (a_box.bigEnd(i) - a_box.smallEnd(i) + 1);
    }

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
    memcpy(data,m_data,m_size*sizeof(double));
#if 0
    for(int i = 0; i < m_size; i++)
    {
        data[i] = m_data[i];
    }
#endif
}

void FArrayBox::copyFromMemory(double *data)
{
    memcpy(m_data,data,m_size*sizeof(double));
#if 0
    for(int i = 0; i < m_size; i++)
    {
        m_data[i] = data[i];
    }
#endif
}

void FArrayBox::set_to_value(double &value)
{
    // Set value.
    for (int i = 0; i < m_size; i++)
    {
        m_data[i] = value;
    }
}

void FArrayBox::set_to_nan()
{
    double snan,big_number;
    set_snan(snan);
    set_big_number(big_number);
    // Set value.
    for (int i = 0; i < m_size; i++)
    {
        m_data[i] = big_number;
    }
}

void FArrayBox::set_to_big_number()
{
    double big_number;
    set_big_number(big_number);
    // Set value.
    for (int i = 0; i < m_size; i++)
    {
        m_data[i] = big_number;
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
Box::Box(const int ll[], const int ur[], const int box_dim):
    m_box_dim(box_dim),
    m_ll(ll,ll+box_dim),
    m_ur(ur,ur+box_dim)
{
}

int Box::smallEnd(int idir) const
{
    return m_ll[idir];
}

int Box::bigEnd(int idir) const
{
    return m_ur[idir];
}

int Box::boxDim() const
{
    return m_box_dim;
}
