#ifndef __AMR_UTILS__
#define __AMR_UTILS__

#define CH_assert(a)

#define CH_SPACEDIM 2

int SpaceDim = 2;

#define Max(a,b)  a

typedef double Real;

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
