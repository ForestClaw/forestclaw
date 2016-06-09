      subroutine geoclaw_average2coarse(mx,my,mbc,meqn,qcoarse,
     &      qfine,maux,aux_coarse,aux_fine,mbathy,
     &      p4est_refineFactor,refratio,igrid)
      implicit none

      integer mx,my,mbc,meqn,p4est_refineFactor,refratio
      integer igrid,maux,mbathy
      double precision qcoarse(meqn,1-mbc:mx+mbc,1-mbc:my+mbc)
      double precision qfine(meqn,1-mbc:mx+mbc,1-mbc:my+mbc)
      double precision aux_coarse(maux,1-mbc:mx+mbc,1-mbc:my+mbc)
      double precision aux_fine(maux,1-mbc:mx+mbc,1-mbc:my+mbc)

      integer i,j, ig, jg, ic_add, jc_add, ii, jj, ifine, jfine
      integer mq
      double precision sum

c     # This should be refratio*refratio.
      integer i1,j1,r2,m
      integer rr2
      parameter(rr2 = 4)
      integer i2(0:rr2-1),j2(0:rr2-1)
      double precision kc, kf, qf

      r2 = refratio*refratio
      if (r2 .ne. rr2) then
         write(6,*) 'average_face_ghost (claw2d_utils.f) ',
     &         '  Refratio**2 is not equal to rr2'
         stop
      endif


c     # Get (ig,jg) for grid from linear (igrid) coordinates
      ig = mod(igrid,refratio)
      jg = (igrid-ig)/refratio

c     # Get rectangle in coarse grid for fine grid.
      ic_add = ig*mx/p4est_refineFactor
      jc_add = jg*my/p4est_refineFactor

      r2 = refratio * refratio

      do mq = 1,meqn
c        # First loop over quadrant (i1,i2)x(j1,j2) of the coarse grid
         do j = 1,my/p4est_refineFactor
            do i = 1,mx/p4est_refineFactor
               i1 = i+ic_add
               j1 = j+jc_add
               m = 0
               do jj = 1,refratio
                  do ii = 1,refratio
                     i2(m) = (i-1)*refratio + ii
                     j2(m) = (j-1)*refratio + jj
                     m = m + 1
                  enddo
               enddo
               if (mq .eq. 1) then
c                 # Average sea surface height rather than just the
c                 # water column height.
                  sum = 0
                  do m = 0,r2-1
                     qf = qfine(mq,i2(m),j2(m)) + 
     &                    aux_fine(mbathy,i2(m),j2(m))
                     sum = sum + qf
                  enddo
                  qcoarse(mq,i1,j1) = sum/r2 - aux_coarse(mq,i1,j1)
               else
c                 # Average momentum components in the usual way.
c                 # But then make sure that no new extrema are created.
                  sum = 0
                  do m = 0,r2-1
                     qf = qfine(mq,i2(m),j2(m))
                     sum = sum + qf
                  enddo
                  qcoarse(mq,i1,j1) = sum/r2
cc                 # check to make sure we are not creating any new extrema
c                  do ii = -1,1
c                     do jj = -1,1
c                        uc(ii,jj) = qcoarse(mq,ic+ii,jc+jj)/
c     &                        qcoarse(1,ic+ii,jc+jj)
c                     enddo
c                  enddo
c
c                  coarseumax = -1d99
c                  coarseumin = 1d99
c                  do ii = -1,1
c                     do jj = -1,1
c                        coarseumax = max(coarseumax,uc(ii,jj))
c                        coarseumin = min(coarseumin,uc(ii,jj))
c                     enddo
c                  enddo
c
c                  redefine = .false.
c                  do ii = 1,refratio
c                     do jj = 1,refratio
c                        iff = (i-1)*refratio + ii
c                        jf = (j-1)*refratio + jj
c                        uf = qfine(mq,iff,jf)/qfine(1,iff,jf)
c                        if (uf .gt. coarseumax .or. uf .lt. coarseumin)
c     &                        then
c                           redefine = .true.
c                        endif
c                     enddo
c                  enddo
c
c                  if (redefine) then
c                     do ii = 1,refratio
c                        do jj = 1,refratio
c                           iff = (i-1)*refratio + ii
c                           jf = (j-1)*refratio + jj
c                           qfine(mq,iff,jf) = qfine(1,iff,jf)*uc(0,0)
c                        enddo
c                     enddo
c                  endif

               endif
            enddo
         enddo
      enddo

      end
