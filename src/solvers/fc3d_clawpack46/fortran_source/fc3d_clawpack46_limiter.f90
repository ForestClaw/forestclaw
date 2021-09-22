subroutine clawpack46_limiter(maxm,meqn,mwaves,mbc,mx, & 
           wave,s,mthlim)

    !! # Apply a limiter to the waves.
    !! # The limiter is computed by comparing the 2-norm of each wave with
    !! # the projection of the wave from the interface to the left or
    !! # right onto the current wave.  For a linear system this would
    !! # correspond to comparing the norms of the two waves.  For a
    !! # nonlinear problem the eigenvectors are not colinear and so the
    !! # projection is needed to provide more limiting in the case where the
    !! # neighboring wave has large norm but points in a different direction
    !! # in phase space.
    !!
    !! # The specific limiter used in each family is determined by the
    !! # value of the corresponding element of the array mthlim, as used in
    !! # the function philim.
    !! # Note that a different limiter may be used in each wave family.
    !!
    !! # dotl and dotr denote the inner product of wave with the wave to
    !! # the left or right.  The norm of the projections onto the wave are then
    !! # given by dotl/wnorm2 and dotr/wnorm2, where wnorm2 is the 2-norm
    !! # of wave.

    implicit none

    integer :: maxm, meqn, mwaves, mbc, mx

    integer :: mthlim(mwaves)
    double precision :: wave(meqn,mwaves,1-mbc:maxm+mbc)
    double precision :: s(mwaves,1-mbc:maxm+mbc)

    integer :: mw,i, m
    double precision :: dotr, wnorm2, dotl, wlimitr
    double precision :: philim

    wave_loop : do mw = 1,mwaves
        if (mthlim(mw) .eq. 0) then
            cycle
        endif
        dotr = 0.d0
        i_loop : do  i = 0, mx+1
            wnorm2 = 0.d0
            dotl = dotr
            dotr = 0.d0
            do m = 1,meqn
                wnorm2 = wnorm2 + wave(m,mw,i)**2
                dotr = dotr + wave(m,mw,i)*wave(m,mw,i+1)
            enddo

            if (i .gt. 0 .and. wnorm2 .ne. 0) then
                if (s(mw,i) .gt. 0.d0) then
                    wlimitr = philim(wnorm2, dotl, mthlim(mw))
                else
                    wlimitr = philim(wnorm2, dotr, mthlim(mw))
                endif

                do m = 1,meqn
                    wave(m,mw,i) = wlimitr * wave(m,mw,i)
                end do
            endif
        end do i_loop
    end do wave_loop

    return
end subroutine  clawpack46_limiter
