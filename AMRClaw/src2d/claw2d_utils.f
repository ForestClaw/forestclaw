      subroutine inputparms(mx_leaf,my_leaf,initial_dt, tfinal,
     &      max_cfl, desired_cfl,nout, src_term, verbose, mcapa,
     &      maux,meqn,mwaves,
     &      maxmwaves,mthlim,mbc,mthbc,order,maxlevel)
      implicit none

      integer mx_leaf, my_leaf
      double precision initial_dt, tfinal, max_cfl, desired_cfl
      integer nout, src_term, mcapa, maux, meqn, mwaves
      integer maxmwaves,mbc, verbose
      integer mthbc(4),mthlim(maxmwaves), order(2)
      integer maxlevel

      integer mw, m


      open(55,file='claw2ez.data')

      read(55,*) mx_leaf
      read(55,*) my_leaf

c     timestepping variables
      read(55,*) nout
      read(55,*) tfinal

      read(55,*) initial_dt
      read(55,*) max_cfl
      read(55,*) desired_cfl

      read(55,*) (order(m),m=1,2)
      read(55,*) verbose
      read(55,*) src_term
      read(55,*) mcapa
      read(55,*) maux

      read(55,*) meqn
      read(55,*) mwaves

      if (mwaves .gt. maxmwaves) then
         write(6,*) 'ERROR : (claw_utils.f) mwaves > maxmwaves'
         write(6,*) 'mwaves = ', mwaves
         write(6,*) 'maxmwaves = ', maxmwaves
         stop
      endif

      read(55,*) (mthlim(mw), mw=1,mwaves)

      read(55,*) mbc
      read(55,*) mthbc(1)
      read(55,*) mthbc(2)
      read(55,*) mthbc(3)
      read(55,*) mthbc(4)

      read(55,*) maxlevel

      close(55)
      end


c     Average fine grid to coarse grid or copy neighboring coarse grid
      subroutine average_ghost_step1(mx,my,mbc,meqn,
     &      qcoarse,qfine,idir,iside,refratio,igrid)
      implicit none

      integer mx,my,mbc,meqn,refratio,igrid,idir,iside
      double precision qcoarse(1-mbc:mx+mbc,1-mbc:my+mbc,meqn)
      double precision qfine(1-mbc:mx+mbc,1-mbc:my+mbc,meqn)
      double precision sum

      integer i,j,ibc,jbc,i1,j1,ii,jj,mq,r2

      r2 = refratio*refratio

      if (refratio .eq. 1) then
c        # We only have to consider the case of exchanges on
c        # the high side of the coarse grid.
c        # We only have to copy from one grid to the other.
         if (idir .eq. 0) then
            do j = 1,my
               do ibc = 1,mbc
                  do mq = 1,meqn
                     qcoarse(mx+ibc,j,mq) = qfine(ibc,j,mq)
                  enddo
               enddo
            enddo
         else
            do i = 1,mx
               do jbc = 1,mbc
                  do mq = 1,meqn
                     qcoarse(i,my+jbc,mq) = qfine(i,jbc,mq)
                  enddo
               enddo
            enddo
         endif
      else
c        # Average fine grid onto coarse grid
         if (idir .eq. 0) then
            do j = (igrid-1)*my/refratio,igrid*my/refratio
               do ibc = 1,mbc
                  do mq = 1,meqn
                     sum = 0
                     do ii = 1,refratio
                        do jj = 1,refratio
                           i1 = (ibc-1)*refratio + ii
                           j1 = (j-1)*refratio + jj
                           if (iside .eq. 0) then
                              sum = sum + qfine(mx-i1+1,j1,mq)
                           else
                              sum = sum + qfine(i1,j1,mq)
                           endif
                        enddo
                     enddo
                     if (iside .eq. 0) then
                        qcoarse(1-ibc,j,mq) = sum/r2
                     else
                        qcoarse(mx+ibc,j,mq) = sum/r2
                     endif
                  enddo
               enddo
            enddo
         else
            do i = (igrid-1)*mx/refratio+1,igrid*mx/refratio
               do jbc = 1,mbc
                  do mq = 1,meqn
                     sum = 0
                     do ii = 1,refratio
                        do jj = 1,refratio
                           i1 = (i-1)*refratio + ii
                           j1 = (jbc-1)*refratio + jj
                           if (iside .eq. 0) then
                              sum = sum + qfine(i1,my-j1+1,mq)
                           else
                              sum = sum + qfine(i1,j1,mq)
                           endif
                        enddo
                     enddo
                     if (iside .eq. 0) then
                        qcoarse(i,1-jbc,mq) = sum/r2
                     else
                        qcoarse(i,my+jbc,mq) = sum/r2
                     endif
                  enddo
               enddo
            enddo
         endif
      endif



      end
