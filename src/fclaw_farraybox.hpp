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

#ifndef FCLAW_FARRAYBOX_HPP
#define FCLAW_FARRAYBOX_HPP

#include <fclaw2d_defs.h>
#include <vector>

void fclaw_farraybox_set_to_nan(double& f);

class Box
{
public:
    Box() = default;
    Box(const int ll[], const int ur[], const int box_dim);
    int smallEnd(int idir) const;
    int bigEnd(int idir) const;
    int boxDim() const;

private:
    int m_box_dim = 0;
    std::vector<int> m_ll;
    std::vector<int> m_ur;
};

class FArrayBox
{
public:

    FArrayBox();
    FArrayBox(const FArrayBox& A);
    ~FArrayBox();
    void define(const Box& a_box, int a_fields);
    double* dataPtr();
    Box box();
    int fields();
    void set_to_value(double &value);
    void set_to_nan();
    void set_to_big_number();
    int size();
    void operator=(const FArrayBox& fbox);
    void copyToMemory(double *data);
    void copyFromMemory(double *data);
private:
    double *m_data;
    int m_size;
    Box m_box;
    int m_fields;
    void set_dataPtr(int size);

};

#endif
