#ifndef __AMR_UTILS__
#define __AMR_UTILS__

#define CH_assert(a)

#define CH_SPACEDIM 2

int SpaceDim = 2;

#define Max(a,b)  a

typedef double Real;


extern "C"
{
    void inputparms_(Real& initial_dt,
                     Real& tfinal,
                     Real& max_cfl,
                     int& nout,
                     int& src_term,
                     int& mcapa,
                     int& maux,
                     int& meqn,
                     int& mwaves,
                     int& maxmwaves,
                     int  *mthlim,
                     int& mbc,
                     int  *mthbc,
                     int  *order);
}

class global_parms
{
public:
    double m_initial_dt;
    double m_tfinal;
    double m_max_cfl;
    int m_nout;
    int m_src_term;
    int m_mcapa;
    int m_maux;
    int m_meqn;
    int m_mwaves;
    int m_maxmwaves;
    int *m_mthlim;
    int m_mbc;
    int *m_mthbc;
    int *m_order;

    global_parms()
        {
            int m_maxmwaves = 20;
            m_mthlim = new int[m_maxmwaves];
            m_mthbc = new int[4];
            m_order = new int[2];
        }

    ~global_parms()
        {
            delete [] m_mthlim;
            delete [] m_mthbc;
            delete [] m_order;
        }

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
