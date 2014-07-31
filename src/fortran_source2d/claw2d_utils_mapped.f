      subroutine average_to_coarse_mapped(mx,my,mbc,meqn,qcoarse,qfine,
     &      areacoarse, areafine, p4est_refineFactor,
     &      refratio, igrid)
      implicit none

      integer mx,my,mbc,meqn,p4est_refineFactor, refratio, igrid
      double precision qcoarse(1-mbc:mx+mbc,1-mbc:my+mbc,meqn)
      double precision   qfine(1-mbc:mx+mbc,1-mbc:my+mbc,meqn)

      double precision areacoarse(-mbc:mx+mbc+1,-mbc:my+mbc+1)
      double precision   areafine(-mbc:mx+mbc+1,-mbc:my+mbc+1)

      integer i,j, ig, jg, ic_add, jc_add, ii, jj, ifine, jfine
      integer mq
      double precision sum, r2, qf, kf, kc

      write(6,*) 'average_to_coarse_mapped (claw2d_utils_mapped) : ',
     &      'This should not get called'
      stop

c     # Get (ig,jg) for grid from linear (igrid) coordinates
      ig = mod(igrid,refratio)
      jg = (igrid-ig)/refratio

c     # Get rectangle in coarse grid for fine grid.
      ic_add = ig*mx/p4est_refineFactor
      jc_add = jg*mx/p4est_refineFactor

c     # Do an area weighted average
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
                     kf = areafine(ifine,jfine)
                     sum = sum + kf*qf
                  enddo
               enddo
               kc = areacoarse(i+ic_add,j+jc_add)
               qcoarse(i+ic_add,j + jc_add,mq) = sum/kc
            enddo
         enddo
      enddo
      end


c     # average ghost cells from 'igrid' neighbor 'qfine' (igrid = 0,1)
c     # to 'qcoarse' at face 'iside'  in direction 'idir' of 'qcoarse'
      subroutine average_face_ghost_mapped(mx,my,mbc,meqn,qcoarse,
     &      qfine,areacoarse, areafine,
     &      idir,iface_coarse,p4est_refineFactor,refratio,igrid)
      implicit none

      integer mx,my,mbc,meqn,refratio,igrid,idir,iface_coarse
      integer p4est_refineFactor
      double precision qfine(1-mbc:mx+mbc,1-mbc:my+mbc,meqn)
      double precision qcoarse(1-mbc:mx+mbc,1-mbc:my+mbc,meqn)
      double precision areacoarse(-mbc:mx+mbc+1,-mbc:my+mbc+1)
      double precision areafine(-mbc:mx+mbc+1,-mbc:my+mbc+1)

      double precision sum, kf, kc, qc,qf

      integer mq,r2
      integer i, ic_add, ibc, ii, ifine
      integer j, jc_add, jbc, jj, jfine

      write(6,*) 'average_to_coarse_mapped (claw2d_utils_mapped) : ',
     &      'This should not get called'
      stop


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
                           kf = areafine(mx-ifine+1,jfine)
                        else
                           qf = qfine(ifine,jfine,mq)
                           kf = areafine(ifine,jfine)
                        endif
                        sum = sum + qf*kf
                     enddo
                  enddo
                  if (iface_coarse .eq. 0) then
                     kc = areacoarse(1-ibc,j+jc_add)
                     qcoarse(1-ibc,j+jc_add,mq) = sum/kc
                  else
                     kc = areacoarse(mx+ibc,j+jc_add)
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
                           kf = areafine(ifine,my-jfine+1)
                        else
                           qf = qfine(ifine,jfine,mq)
                           kf = areafine(ifine,jfine)
                        endif
                        sum = sum + kf*qf
                     enddo
                  enddo
                  if (iface_coarse .eq. 2) then
                     kc = areacoarse(i+ic_add,1-jbc)
                     qcoarse(i+ic_add,1-jbc,mq) = sum/kc
                  else
                     kc = areacoarse(i+ic_add,my+jbc)
                     qcoarse(i+ic_add,my+jbc,mq) = sum/kc
                  endif
               enddo
            enddo
         endif
      enddo

      end

      subroutine fixcapaq2(mx,my,mbc,meqn,qcoarse,qfine,
     &      areacoarse,areafine, p4est_refineFactor,
     &      refratio,igrid)
      implicit none

      integer mx,my,mbc,meqn, refratio, igrid
      integer p4est_refineFactor

      double precision qcoarse(1-mbc:mx+mbc,1-mbc:my+mbc,meqn)
      double precision qfine(1-mbc:mx+mbc,1-mbc:my+mbc,meqn)
      double precision areacoarse(-mbc:mx+mbc+1,-mbc:my+mbc+1)
      double precision   areafine(-mbc:mx+mbc+1,-mbc:my+mbc+1)

      integer i,j,ii, jj, ifine, jfine, m, ig, jg, ic_add, jc_add
      double precision kf, kc, r2, sum, cons_diff, qf, qc

      write(6,*) 'average_to_coarse_mapped (claw2d_utils_mapped) : ',
     &      'This should not get called'
      stop


c     # Get (ig,jg) for grid from linear (igrid) coordinates
      ig = mod(igrid,refratio)
      jg = (igrid-ig)/refratio

c     # Get rectangle in coarse grid for fine grid.
      ic_add = ig*mx/p4est_refineFactor
      jc_add = jg*my/p4est_refineFactor

c     # ------------------------------------------------------
c     # This routine ensures that the interpolated solution
c     # has the same mass as the coarse grid solution
c     # -------------------------------------------------------


      r2 = refratio*refratio
      do m = 1,meqn
         do i = 1,mx/p4est_refineFactor
            do j = 1,my/p4est_refineFactor
               sum = 0.d0
               do ii = 1,refratio
                  do jj = 1,refratio
                     ifine = (i-1)*refratio + ii
                     jfine = (j-1)*refratio + jj
                     kf = areafine(ifine,jfine)
                     qf = qfine(ifine,jfine,m)
                     sum = sum + kf*qf
                  enddo
               enddo

               kc = areacoarse(i+ic_add,j+jc_add)
               qc = qcoarse(i+ic_add, j+jc_add,m)
               cons_diff = (qc*kc - sum)/r2

               do ii = 1,refratio
                  do jj = 1,refratio
                     ifine  = (i-1)*refratio + ii
                     jfine  = (j-1)*refratio + jj
                     kf = areafine(ifine,jfine)
                     if (kf .eq. 0) then
                        write(6,*) 'fixcapaq : kf = 0'
                        stop
                     endif
                     qfine(ifine,jfine,m) = qfine(ifine,jfine,m) +
     &                     cons_diff/kf
                  enddo
               enddo
            enddo  !! end of meqn
         enddo
      enddo

      end
