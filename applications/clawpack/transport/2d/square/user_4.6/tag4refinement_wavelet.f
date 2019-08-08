      subroutine square_tag4refinement_wavelet(mx,my,mbc,
     &      meqn, xlower,ylower,dx,dy,blockno,
     &      q, tag_threshold, init_flag,tag_patch)
      implicit none

      integer mx,my, mbc, meqn, tag_patch, init_flag
      integer blockno
      double precision xlower, ylower, dx, dy
      double precision tag_threshold
      double precision q(1-mbc:mx+mbc,1-mbc:my+mbc,meqn)

      integer initchoice
      common /initchoice_comm/ initchoice

      integer*8 cont, get_context
      logical fclaw2d_map_is_used

      integer i,j, mq
      double precision qmin, qmax, xc, yc
      double precision xp,yp,zp

      double precision qcoarse(1-mbc:mx+mbc, 1-mbc:my+mbc,meqn)
      double precision qfine(1-mbc:mx+mbc, 1-mbc:my+mbc,meqn)
      double precision areacoarse(-mbc:mx+mbc+1,-mbc:my+mbc+1)
      double precision   areafine(-mbc:mx+mbc+1,-mbc:my+mbc+1)
      integer manifold
      
      tag_patch = 0

      cont = get_context()

      manifold = 0
      call wavelet_average2coarse(mx,my,mbc,meqn,
     &      q,qcoarse, areacoarse, areafine, manifold)

      call wavelet_interpolate2fine3
     &     (mx,my,mbc,meqn,qcoarse, qfine, areacoarse, 
     &      areafine, manifold)


c     # Refine based only on first variable in system.
      mq = 1
      qmin = q(1,1,mq)
      qmax = q(1,1,mq)
      do j = 1,my
         do i = 1,mx
             if (abs(q(i,j,1) - qfine(i,j,1)) .gt. tag_threshold) then
                 tag_patch = 1
                 return
            endif
         enddo
      enddo

      end



      subroutine wavelet_average2coarse(mx,my,mbc,meqn,
     &      qfine,qcoarse,areacoarse, areafine, manifold)
      implicit none

      integer mx,my,mbc,meqn
      integer manifold
      double precision qcoarse(1-mbc:mx+mbc,1-mbc:my+mbc,meqn)
      double precision qfine(1-mbc:mx+mbc,1-mbc:my+mbc,meqn)

c     # these will be empty if we are not on a manifold.
      double precision areacoarse(-mbc:mx+mbc+1,-mbc:my+mbc+1)
      double precision   areafine(-mbc:mx+mbc+1,-mbc:my+mbc+1)

      integer i,j, mq
      double precision sum
      logical is_manifold

c     # This should be refratio*refratio.
      integer ii,jj,m
      integer i2(0:3),j2(0:3)
      double precision kc, kf, qf

      is_manifold = manifold .eq. 1

c     # 'iface' is relative to the coarse grid

      do mq = 1,meqn
         do j = 0,my/2+1
            do i = 0,mx/2+1
               m = 0
               do jj = 1,2
                  do ii = 1,2
                     i2(m) = (i-1)*2 + ii
                     j2(m) = (j-1)*2 + jj
                     m = m + 1
                  enddo
               enddo
               if (is_manifold) then
                  sum = 0
                  do m = 0,3
                     qf = qfine(i2(m),j2(m),mq)
                     kf = areafine(i2(m),j2(m))
                     sum = sum + kf*qf
                  enddo
                  kc = areacoarse(i,j)
                  qcoarse(i,j,mq) = sum/kc
               else
                  sum = 0
                  do m = 0,3
                     qf = qfine(i2(m),j2(m),mq)
                     sum = sum + qf
                  enddo
                  qcoarse(i,j,mq) = sum/4
               endif
            enddo
         enddo
      enddo
      end


      subroutine wavelet_interpolate2fine3
     &     (mx,my,mbc,meqn,qcoarse, qfine, areacoarse, 
     &      areafine, manifold)
      implicit none

      integer mx,my,mbc,meqn
      integer manifold

      double precision qcoarse(1-mbc:mx+mbc,1-mbc:my+mbc,meqn)
      double precision qfine(1-mbc:mx+mbc,1-mbc:my+mbc,meqn)

      double precision areacoarse(-mbc:mx+mbc+1,-mbc:my+mbc+1)
      double precision   areafine(-mbc:mx+mbc+1,-mbc:my+mbc+1)

      integer i,j,ii,jj,ic,jc, mq, mth
      double precision qc, shiftx, shifty, sl, sr, gradx, grady
      double precision fclaw2d_clawpatch_compute_slopes
      integer order

      double precision qv(8), value, pv(5), RinvQt(5,8)
      integer iff,jff, k

      data RinvQt /  -1.6666666666666663d-01, 1.6666666666666663d-01,
     &      0.1d0, -0.25d0, 0.1d0, 0, 1.6666666666666663d-01, 
     &      -0.2d0, 0, 0.3d0, 1.6666666666666663d-01, 
     &      1.6666666666666663d-01,
     &      0.1d0, 0.25d0, 0.1d0, -1.6666666666666663d-01,
     &      0, 0.3d0, 0, -0.2d0, 1.6666666666666663d-01,
     &      0, 0.3d0, 0, -0.2d0, -1.6666666666666663d-01,
     &     -1.6666666666666663d-01, 0.1d0, 0.25d0,
     &      0.1d0, 0, -1.6666666666666663d-01, -0.2d0,
     &      0, 0.3d0, 1.6666666666666663d-01,
     &     -1.6666666666666663d-01, 0.1d0, -0.25d0, 0.1d0/



c     # Use limiting done in AMRClaw.
      mth = 5

      order = 3
c     # Only averaged to the lower corner of the coarse grid
      do mq = 1,meqn
         do j = 1,my/2
            do i = 1,mx/2
               ic = i 
               jc = j
               qc = qcoarse(ic,jc,mq)

               if (order .eq. 2) then
c                  # Compute limited slopes in both x and y. Note we are not
c                  # really computing slopes, but rather just differences.
c                  # Scaling is accounted for in 'shiftx' and 'shifty', below.
                   sl = (qc - qcoarse(ic-1,jc,mq))
                   sr = (qcoarse(ic+1,jc,mq) - qc)
                   gradx = fclaw2d_clawpatch_compute_slopes(sl,sr,mth)

                   sl = (qc - qcoarse(ic,jc-1,mq))
                   sr = (qcoarse(ic,jc+1,mq) - qc)
                   grady = fclaw2d_clawpatch_compute_slopes(sl,sr,mth)

                   do ii = 1,2
                       do jj = 1,2
                           shiftx = (ii - 1 - 0.5d0)/2
                           shifty = (jj - 1 - 0.5d0)/2
                           qfine((i-1)*2 + ii,(j-1)*2 + jj,mq)
     &                           = qc + shiftx*gradx + shifty*grady
                       end do
                  end do
               elseif (order .eq. 3) then

c                  # This is the right hand side of the LLSQ problem
                   qv(1) = qcoarse(ic-1,jc+1,mq) - qc
                   qv(2) = qcoarse(ic,  jc+1,mq) - qc 
                   qv(3) = qcoarse(ic+1,jc+1,mq) - qc 
                   qv(4) = qcoarse(ic-1,jc  ,mq) - qc 
                   qv(5) = qcoarse(ic+1,jc  ,mq) - qc 
                   qv(6) = qcoarse(ic-1,jc-1,mq) - qc 
                   qv(7) = qcoarse(ic,  jc-1,mq) - qc 
                   qv(8) = qcoarse(ic+1,jc-1,mq) - qc 
      
                   do ii = 1,5
                       pv(ii) = 0
                       do k = 1,8
                           pv(ii) = pv(ii) + RinvQt(ii,k)*qv(k)
                       end do
                   end do


c                  # Fill in refined values on coarse grid cell (ic,jc)
                   do ii = 1,2
                       do jj = 1,2
                           shiftx = (ii - 1 - 0.5d0)/2
                           shifty = (jj - 1 - 0.5d0)/2
                           iff = (i-1)*2 + ii
                           jff = (j-1)*2 + jj
                           value = pv(1)*shiftx + pv(2)*shifty + 
     &                            pv(4)*shiftx*shifty
                           qfine(iff,jff,mq) = qc + value
                       enddo
                   enddo
               endif
            end do
         enddo
      enddo

      end







