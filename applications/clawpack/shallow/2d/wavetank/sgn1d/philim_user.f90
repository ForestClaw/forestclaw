
double precision function philim_user(a,b,meth)
    implicit none

    double precision a,b
    integer meth

    double precision wlimitr, th, r

!!     # Compute a limiter based on wave strengths a and b.
!!     # meth determines what limiter is used.
!!     # a is assumed to be nonzero.

!!     # NOTE: This routine is obsolete.  Instead of using limiter.f,
!!     # which calls philim.f for every wave, it is more efficient to
!!     # use inlinelimiter.f, which eliminates all these function calls
!!     # to philim.  If you wish to change the limiter function and are
!!     # using inlinelimiter.f, the formulas must be changed in that routine.

    r = b/a

    if (meth .eq. 6) then
        !! ------------------------------
        !! # Sweby
        !! ------------------------------
        wlimitr = dmax1(0.d0, dmin1(r,1.5d0),dmin1(1.5d0*r,1.d0))

    else if (meth .eq. 7) then
        !! ------------------------------
        !! # Generalized minmod
        !! ------------------------------
        th = 1.3
        wlimitr = dmax1(0.d0, dmin1(th*r,(1+r)/2.d0,th));
    endif

    philim_user = wlimitr

    return
end function philim_user
