      subroutine inputparms(mx_leaf,my_leaf,initial_dt, tfinal,
     &      max_cfl, desired_cfl,nout, src_term, verbose, mcapa,
     &      maux,meqn,mwaves,
     &      maxmwaves,mthlim,mbc,mthbc,order,minlevel,maxlevel,
     &      icycle)
      implicit none

      integer mx_leaf, my_leaf
      double precision initial_dt, tfinal, max_cfl, desired_cfl
      integer nout, src_term, mcapa, maux, meqn, mwaves
      integer maxmwaves,mbc, verbose
      integer mthbc(4),mthlim(maxmwaves), order(2)
      integer maxlevel, minlevel, icycle
      logical subcycle

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

      read(55,*) minlevel
      read(55,*) maxlevel

      read(55,*) subcycle
      if (subcycle) then
         icycle = 1
      else
         icycle = 0
      endif


      close(55)
      end


c     # Exchange edge ghost data with neighboring grid at same level.
      subroutine exchange_face_ghost(mx,my,mbc,meqn, qthis,qneighbor,
     &      idir)
      implicit none

      integer mx,my,mbc,meqn,igrid,idir
      double precision qthis(1-mbc:mx+mbc,1-mbc:my+mbc,meqn)
      double precision qneighbor(1-mbc:mx+mbc,1-mbc:my+mbc,meqn)

      integer i,j,ibc,jbc,mq

c     # We only have to consider the case of exchanges on
c     # the high side of 'this' grid.
c     # We do need to do a complete exchange, though.
      if (idir .eq. 0) then
         do j = 1,my
            do ibc = 1,mbc
               do mq = 1,meqn
c                 # Exchange at high side of 'this' grid in
c                 # x-direction (idir == 0)
                  qthis(mx+ibc,j,mq) = qneighbor(ibc,j,mq)
                  qneighbor(ibc-mbc,j,mq) = qthis(mx-mbc+ibc,j,mq)
               enddo
            enddo
         enddo
      else
         do i = 1,mx
            do jbc = 1,mbc
               do mq = 1,meqn
c                 # Exchange at high side of 'this' grid in
c                 # y-direction (idir == 1)
                  qthis(i,my+jbc,mq) = qneighbor(i,jbc,mq)
                  qneighbor(i,jbc-mbc,mq) = qthis(i,my-mbc+jbc,mq)
               enddo
            enddo
         enddo
      endif
      end



c     # Average fine grid to a coarse grid neighbor or copy from neighboring
c     # grid at same level.
      subroutine ghost_cell_average(mx,my,mbc,meqn,
     &      qfine,qcoarse,idir,iside,refratio,igrid)
      implicit none

      integer mx,my,mbc,meqn,refratio,igrid,idir,iside
      double precision qfine(1-mbc:mx+mbc,1-mbc:my+mbc,meqn)
      double precision qcoarse(1-mbc:mx+mbc,1-mbc:my+mbc,meqn)
      double precision sum

      integer i,j,ibc,jbc,i1,j1,ii,jj,mq,r2

      r2 = refratio*refratio

c     # Average fine grid onto coarse grid
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
                        if (iside .eq. 2) then
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

      end

      subroutine interpolate_face_ghost(mx,my,mbc,meqn,qfine,qcoarse,
     &      idir,iside,refratio,igrid)
      implicit none
      integer mx,my,mbc,meqn,refratio,igrid,idir,iside
      double precision qfine(1-mbc:mx+mbc,1-mbc:my+mbc,meqn)
      double precision qcoarse(1-mbc:mx+mbc,1-mbc:my+mbc,meqn)
      double precision sum

      integer i,j,ibc,jbc,i1,j1,ii,jj,mq,r2


      end


c     Average fine grid to coarse grid or copy neighboring coarse grid
      subroutine average_corner(mx,my,mbc,meqn,
     &      qcoarse,qfine,icorner,refratio)
      implicit none

      integer mx,my,mbc,meqn,refratio,icorner
      double precision qcoarse(1-mbc:mx+mbc,1-mbc:my+mbc,meqn)
      double precision qfine(1-mbc:mx+mbc,1-mbc:my+mbc,meqn)
      double precision sum

      integer i,j,ibc,jbc,i1,j1,ii,jj,mq,r2

      r2 = refratio*refratio

      if (refratio .eq. 1) then
c        # We only need to worry about corners 1 and 3 (lr and ur).
c        # for the  complete exchange.  The other corners are somebody
c        # else's (lr,ur) corners.
         do mq = 1,meqn
            do ibc = 1,mbc
               do jbc = 1,mbc
c                 # Exchange corner information at boundaries.
                  if (icorner .eq. 1) then
c                    # Fix this!
                     qcoarse(mx+ibc,jbc-mbc,mq) =
     &                     qfine(ibc,my+jbc-mbc,mq)
                     qfine(ibc-mbc,my+jbc,mq) =
     &                     qcoarse(mx+ibc-mbc,jbc,mq)
                  elseif (icorner .eq. 3) then
                     qcoarse(mx+ibc,my+jbc,mq) =
     &                     qfine(ibc,jbc,mq)
                     qfine(ibc-mbc,jbc-mbc,mq) =
     &                     qcoarse(mx+ibc-mbc,my+jbc-mbc,mq)
                  endif
               enddo
            enddo
         enddo
      else
c        # Average fine grid onto coarse grid
         write(6,'(A,A)') 'average_corner_step1 : fine grid ',
     &         ' averaging at corners not yet implemented'
         stop
         do mq = 1,meqn
            do ibc = 1,mbc
               do jbc = 1,mbc
                  if (icorner .eq. 0) then
                  elseif (icorner .eq. 1) then
                  elseif (icorner .eq. 2) then
                  else
                  endif
               enddo
            enddo
         enddo
      endif

      end


      subroutine set_phys_corner_ghost(mx,my,mbc,meqn,q,icorner,t,dt,
     &      mthbc)
      implicit none

      integer mx,my,mbc,meqn,icorner, mthbc(4)
      double precision t,dt
      double precision q(1-mbc:mx+mbc,1-mbc:my+mbc,meqn)
      integer i,j

c     # Do something here....

      end


      subroutine exchange_phys_corner_ghost(mx,my,mbc,meqn,
     &      qthis, qneighbor, icorner, iside)
      implicit none

      integer mx, my, mbc, meqn, iside, icorner
      double precision qthis(1-mbc:mx+mbc,1-mbc:my+mbc,meqn)
      double precision qneighbor(1-mbc:mx+mbc,1-mbc:my+mbc,meqn)

      integer ibc, jbc, mq

c     # Fill in corner ghost cells that overlap the physical boundary. In this
c     case, the corner ghost cells are copied from a face neighbor.
      do mq = 1,meqn
         do ibc = 1,mbc
            do jbc = 1,mbc
               if (iside .eq. 1) then
                  if (icorner .eq. 1) then
                     qthis(mx+ibc,jbc-mbc,mq) =
     &                     qneighbor(ibc,jbc-mbc,mq)
                     qneighbor(ibc-mbc,jbc-mbc,mq) =
     &                     qthis(mx+ibc-mbc,jbc-mbc,mq)
                  elseif(icorner .eq. 3) then
                     qthis(mx+ibc,my+jbc,mq) =
     &                     qneighbor(ibc,my+jbc,mq)
                     qneighbor(ibc-mbc,my+jbc,mq) =
     &                     qthis(mx+ibc-mbc,my+jbc,mq)
                  endif
               elseif (iside .eq. 3) then
                  if (icorner .eq. 2) then
                     qthis(ibc-mbc,my+jbc,mq) =
     &                     qneighbor(ibc-mbc,jbc,mq)
                     qneighbor(ibc-mbc,jbc-mbc,mq) =
     &                     qthis(ibc-mbc,my+jbc-mbc,mq)
                  elseif(icorner .eq. 3) then
                     qthis(mx+ibc,my+jbc,mq) =
     &                     qneighbor(mx+ibc,jbc,mq)
                     qneighbor(mx+ibc,jbc-mbc,mq) =
     &                     qthis(mx+ibc,my+jbc-mbc,mq)
                  endif
               endif
            enddo
         enddo
      enddo



      end


      subroutine exchange_corner_ghost(mx,my,mbc,meqn,
     &      qthis, qneighbor, icorner)
      implicit none

      integer mx, my, mbc, meqn, icorner
      double precision qthis(1-mbc:mx+mbc,1-mbc:my+mbc,meqn)
      double precision qneighbor(1-mbc:mx+mbc,1-mbc:my+mbc,meqn)

      integer mq, ibc, jbc

c     # NOTE : qneighbor is not yet a valid pointer.

c     # Only exchanging high side corners

c     # We only need to worry about corners 1 and 3 (lr and ur).
c     # for the  complete exchange.  The other corners are somebody
c     # else's (lr,ur) corners.
      do mq = 1,meqn
         do ibc = 1,mbc
            do jbc = 1,mbc
c              # Exchange corner information at boundaries.
c              # Do this until we get the corner ids fixed
               if (icorner .eq. 1) then
                  qthis(mx+ibc,jbc-mbc,mq) =
     &                  qneighbor(ibc,my+jbc-mbc,mq)
                  qneighbor(ibc-mbc,my+jbc,mq) =
     &                  qthis(mx+ibc-mbc,jbc,mq)
               elseif (icorner .eq. 3) then
                  qthis(mx+ibc,my+jbc,mq) =
     &                  qneighbor(ibc,jbc,mq)
                  qneighbor(ibc-mbc,jbc-mbc,mq) =
     &                  qthis(mx+ibc-mbc,my+jbc-mbc,mq)
               endif
            enddo
         enddo
      enddo


      end
