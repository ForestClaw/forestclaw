      subroutine average_to_coarse_mapped(mx,my,mbc,meqn,qcoarse,qfine,
     &      auxcoarse, auxfine, maux, p4est_refineFactor,
     &      refratio, igrid)
      implicit none

      integer mx,my,mbc,meqn,p4est_refineFactor, refratio, igrid, maux
      double precision qcoarse(1-mbc:mx+mbc,1-mbc:my+mbc,meqn)
      double precision qfine(1-mbc:mx+mbc,1-mbc:my+mbc,meqn)
      double precision auxcoarse(1-mbc:mx+mbc,1-mbc:my+mbc,maux)
      double precision auxfine(1-mbc:mx+mbc,1-mbc:my+mbc,maux)

      integer i,j, ig, jg, ic_add, jc_add, ii, jj, ifine, jfine
      integer mq
      double precision sum, r2, qf, kf, kc

c     # Get (ig,jg) for grid from linear (igrid) coordinates
      ig = mod(igrid,refratio)
      jg = (igrid-ig)/refratio

c     # Get rectangle in coarse grid for fine grid.
      ic_add = ig*mx/p4est_refineFactor
      jc_add = jg*mx/p4est_refineFactor

      r2 = refratio*refratio
      do mq = 1,meqn
         do j = 1,my/p4est_refineFactor
            do i = 1,mx/p4est_refineFactor
               sum = 0
               do ii = 1,refratio
                  do jj = 1,refratio
                     ifine = (i-1)*refratio + ii
                     jfine = (j-1)*refratio + jj
                     qf = qfine(ifine,jfine,mq)
                     kf = auxfine(ifine,jfine,1)
                     sum = sum + kf*qf
                  enddo
               enddo
               kc = r2*auxcoarse(i+ic_add,j+jc_add,1)
               qcoarse(i+ic_add,j + jc_add,mq) = sum/kc
            enddo
         enddo
      enddo
      end


c     # average ghost cells from 'igrid' neighbor 'qfine' (igrid = 0,1)
c     # to 'qcoarse' at face 'iside'  in direction 'idir' of 'qcoarse'
      subroutine average_face_ghost_mapped(mx,my,mbc,meqn,qcoarse,
     &      qfine,auxcoarse, auxfine, maux,
     &      idir,iface_coarse,p4est_refineFactor,refratio,igrid)
      implicit none

      integer mx,my,mbc,meqn,refratio,igrid,idir,iface_coarse
      integer p4est_refineFactor, maux
      double precision qfine(1-mbc:mx+mbc,1-mbc:my+mbc,meqn)
      double precision qcoarse(1-mbc:mx+mbc,1-mbc:my+mbc,meqn)
      double precision auxcoarse(1-mbc:mx+mbc,1-mbc:my+mbc,maux)
      double precision auxfine(1-mbc:mx+mbc,1-mbc:my+mbc,maux)

      double precision sum, kf, kc, qf

      integer mq,r2
      integer i, ic_add, ibc, ii, ifine
      integer j, jc_add, jbc, jj, jfine

c     # 'iface' is relative to the coarse grid

      r2 = refratio*refratio

c     # Average fine grid onto coarse grid
      do mq = 1,meqn
         if (idir .eq. 0) then
            jc_add = igrid*my/p4est_refineFactor
            do j = 1,my/p4est_refineFactor
               do ibc = 1,mbc
c                 # ibc = 1 corresponds to first layer of ghost cells, and
c                 # ibc = 2 corresponds to the second layer
                  sum = 0
                  do ii = 1,refratio
                     do jj = 1,refratio
                        ifine = (ibc-1)*refratio + ii
                        jfine = (j-1)*refratio + jj
                        if (iface_coarse .eq. 0) then
                           qf = qfine(mx-ifine+1,jfine,mq)
                           kf = auxfine(mx-ifine+1,jfine,1)
                        else
                           qf = qfine(ifine,jfine,mq)
                           kf = auxfine(ifine,jfine,1)
                        endif
                        sum = sum + qf*kf
                     enddo
                  enddo
                  if (iface_coarse .eq. 0) then
                     kc = r2*auxcoarse(1-ibc,j+jc_add,1)
                     qcoarse(1-ibc,j+jc_add,mq) = sum/kc
                  else
                     kc = r2*auxcoarse(mx+ibc,j+jc_add,1)
                     qcoarse(mx+ibc,j+jc_add,mq) = sum/kc
                  endif
               enddo
            enddo
         else
            ic_add = igrid*mx/p4est_refineFactor
            do i = 1,mx/p4est_refineFactor
               do jbc = 1,mbc
                  sum = 0
                  do ii = 1,refratio
                     do jj = 1,refratio
                        ifine = (i-1)*refratio + ii
                        jfine = (jbc-1)*refratio + jj
                        if (iface_coarse .eq. 2) then
                           qf = qfine(ifine,my-jfine+1,mq)
                           kf = auxfine(ifine,my-jfine+1,1)
                        else
                           qf = qfine(ifine,jfine,mq)
                           kf = auxfine(ifine,jfine,1)
                        endif
                        sum = sum + kf*qf
                     enddo
                  enddo
                  if (iface_coarse .eq. 2) then
                     kc = r2*auxcoarse(i+ic_add,1-jbc,1)
                     qcoarse(i+ic_add,1-jbc,mq) = sum/kc
                  else
                     kc = r2*auxcoarse(i+ic_add,my+jbc,1)
                     qcoarse(i+ic_add,my+jbc,mq) = sum/kc
                  endif
               enddo
            enddo
         endif
      enddo

      end

      subroutine fixcapaq2(mx,my,mbc,meqn,qcoarse,qfine,
     &      auxcoarse,auxfine,maux, p4est_refineFactor,
     &      refratio,igrid)
      implicit none

      integer mx,my,mbc,meqn, maux, refratio, igrid
      integer p4est_refineFactor

      double precision qcoarse(1-mbc:mx+mbc,1-mbc:my+mbc,meqn)
      double precision qfine(1-mbc:mx+mbc,1-mbc:my+mbc,meqn)
      double precision auxcoarse(1-mbc:mx+mbc,1-mbc:my+mbc,maux)
      double precision auxfine(1-mbc:mx+mbc,1-mbc:my+mbc,maux)

      integer i,j,ii, jj, ifine, jfine, m, ig, jg, ic_add, jc_add
      double precision kf, kc, r2, sum, cons_diff, qf, qc

c     # Get (ig,jg) for grid from linear (igrid) coordinates
      ig = mod(igrid,refratio)
      jg = (igrid-ig)/refratio

c     # Get rectangle in coarse grid for fine grid.
      ic_add = ig*mx/p4est_refineFactor
      jc_add = jg*my/p4est_refineFactor


      r2 = refratio*refratio
      do m = 1,meqn
         do i = 1,mx/p4est_refineFactor
            do j = 1,my/p4est_refineFactor
               sum = 0.d0
               do ii = 1,refratio
                  do jj = 1,refratio
                     ifine = (i-1)*refratio + ii
                     jfine = (j-1)*refratio + jj
                     kf = auxfine(ifine,jfine,1)
                     sum = sum + kf*qfine(ifine,jfine,m)
                  enddo
               enddo

               kc = auxcoarse(i,j,1)
               qc = qcoarse(i+ic_add, j+jc_add,m)
               cons_diff = qf*kc - sum/r2

               do ii = 1,refratio
                  do jj = 1,refratio
                     ifine  = (i-1)*refratio + ii
                     jfine  = (j-1)*refratio + jj
                     kf = auxfine(ifine,jfine,1)
                     qfine(ifine,jfine,m) = qfine(ifine,jfine,m) +
     &                     cons_diff/kf
                  enddo
               enddo
            enddo  !! end of meqn
         enddo
      enddo

c      write(6,'(A,E24.16)') 'fixcapaq : maxdiff = ',maxdiff
      end


c      subroutine coarsen(mx,my,mbc,meqn,refratio, q, qcoarse)
c      implicit none
c
c      integer mx,my,mbc,meqn, refratio
c      double precision q(1-mbc:mx+mbc,1-mbc:my+mbc,meqn)
c      double precision qcoarse(1-mbc:mx/refratio+mbc,
c     &      1-mbc:my/refratio+mbc,meqn)
c      double precision aux(1-mbc:mx/refratio+mbc,
c     &      1-mbc:my/refratio+mbc,meqn)
c      double precision auxcoarse(1-mbc:mx/refratio+mbc,
c     &      1-mbc:my/refratio+mbc,meqn)
c
c      integer i,j,ii,jj, ifine,jfine, r2, mq
c      double precision sum
c
c      r2 = refratio**2
c
c      do mq = 1,meqn
c         do i = 1,mx/refratio
c            do j = 1,my/refratio
c               sum = 0
c               do ii = 1,refratio
c                  do jj = 1,refratio
c                     ifine = (i - 1)*refratio + ii
c                     jfine = (j - 1)*refratio + jj
c                     sum = sum + q(ifine,jfine,mq)
c                  enddo
c               enddo
c               qcoarse(i,j,mq) = sum/r2
c            enddo
c         enddo
cc     # Coarsen ghost cells
c
c      enddo
c
c
c
c      end
