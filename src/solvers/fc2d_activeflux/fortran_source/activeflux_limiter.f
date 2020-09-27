c
c
c     =====================================================
      subroutine activeflux_limiter(maxm,meqn,mwaves,mbc,mx,
     &      wave,s,mthlim)
c     =====================================================
c
c     # Apply a limiter to the waves.
c     # The limiter is computed by comparing the 2-norm of each wave with
c     # the projection of the wave from the interface to the left or
c     # right onto the current wave.  For a linear system this would
c     # correspond to comparing the norms of the two waves.  For a
c     # nonlinear problem the eigenvectors are not colinear and so the
c     # projection is needed to provide more limiting in the case where the
c     # neighboring wave has large norm but points in a different direction
c     # in phase space.
c
c     # The specific limiter used in each family is determined by the
c     # value of the corresponding element of the array mthlim, as used in
c     # the function philim.
c     # Note that a different limiter may be used in each wave family.
c
c     # dotl and dotr denote the inner product of wave with the wave to
c     # the left or right.  The norm of the projections onto the wave are then
c     # given by dotl/wnorm2 and dotr/wnorm2, where wnorm2 is the 2-norm
c     # of wave.
c
      implicit none

      integer maxm, meqn, mwaves, mbc, mx

      integer  mthlim(mwaves)
      double precision wave(1-mbc:maxm+mbc, meqn, mwaves)
      double precision s(1-mbc:maxm+mbc, mwaves)

      integer mw,i, m, ibc
      double precision dotr, wnorm2, dotl, wlimitr
      double precision philim

c      write(6,*) 'Messed with limiter.f'

c     Memory errors are all coming from boundary conditions not be
c     set correctly in the sphere case.
      do m = 1,meqn
         do mw = 1,mwaves
            do ibc = 1,mbc
c               wave(1-ibc,m,mw) = 0
c               wave(mx+ibc,m,mw) = 0
            enddo
         enddo
      enddo


c
      do mw = 1,mwaves
         if (mthlim(mw) .eq. 0) then
            continue
         endif
         dotr = 0.d0
         do  i = 0, mx+1
            wnorm2 = 0.d0
            dotl = dotr
            dotr = 0.d0
            do m = 1,meqn
               wnorm2 = wnorm2 + wave(i,m,mw)**2
               dotr = dotr + wave(i,m,mw)*wave(i+1,m,mw)
            enddo

            if (i .gt. 0 .and. wnorm2 .ne. 0) then
               if (s(i,mw) .gt. 0.d0) then
                  wlimitr = philim(wnorm2, dotl, mthlim(mw))
               else
                  wlimitr = philim(wnorm2, dotr, mthlim(mw))
               endif

               do m=1,meqn
                  wave(i,m,mw) = wlimitr * wave(i,m,mw)
               enddo
            endif
         enddo
      enddo

      return
      end
