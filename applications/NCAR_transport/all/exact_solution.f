c     # ------------------------------------------------------
c     # EXACT_SOLUTION
c     #
c     # Compute the exact solution, to be used either as an
c     # initial condition (t=0) or for doing convergence tests
c     # later.
c     #
c     # This routine computes the exact solution on the whole
c     # array of values.
c     #
c     # This file assumes that the Gaussians are being advected
c     # along straight trajectories in a periodic fashion.
c     # ------------------------------------------------------
      subroutine exact_solution(mx,my,mbc,xp,yp,zp,q0,qexact,t)
      implicit none

      integer mx,my,mbc
      double precision t
      double precision xp(-mbc:mx+mbc+1,-mbc:my+mbc+1)
      double precision yp(-mbc:mx+mbc+1,-mbc:my+mbc+1)
      double precision zp(-mbc:mx+mbc+1,-mbc:my+mbc+1)
      double precision qexact(0:mx+1,0:my+1)
      double precision q0(0:mx+1,0:my+1)

      double precision v(3), x,y,z
      double precision gaussian_sum
      double precision cosine_bell_sum
      double precision slotted_disk_sum
      double precision kappa, tfinal

      integer i,j, ichoice, get_init_choice

      ichoice = get_init_choice()

      call set_gaussian_locations(t)

      call get_wind_parms(kappa,tfinal)


      do i = 0,mx+1
         do j = 0,my+1
            x = xp(i,j)
            y = yp(i,j)
            z = zp(i,j)
            if (abs(t - tfinal) .lt. 1e-8) then
               qexact(i,j) = q0(i,j)
            else
               if (ichoice .eq. 1) then
                  qexact(i,j) = gaussian_sum(x,y,z)
               elseif (ichoice .eq. 2 .or. ichoice .eq. 3) then
                  qexact(i,j) = cosine_bell_sum(x,y,z)
               elseif (ichoice .eq. 4) then
                  qexact(i,j) = slotted_disk_sum(x,y,z)
               endif
            endif
         enddo
      enddo


      end


      subroutine set_gaussian_locations(t)
      implicit none

      double precision t

      integer i,m
      double precision wc_new(3,2), wr(3)
      integer get_n_locations, n

      n = get_n_locations()

      do i = 1,n
         call get_initial_gaussian_locations(i,wr)
         call rotate_initcond(wr,t)
         do m = 1,3
            wc_new(m,i) = wr(m)
         enddo
      enddo

      call set_td_gaussian_locations(wc_new,n)

      end


      subroutine rotate_initcond(wr,t)
      implicit none

      double precision wr(3),t, th, kappa, tfinal
      double precision vrot(3), r(3,3), pi, uth
      integer i,j,k

      common /compi/ pi

c      th = t*uth
      call get_wind_parms(kappa,tfinal)

      uth = 2*pi/tfinal
      th = t*uth

      do i = 1,3
         do j = 1,3
            r(i,j) = 0
         enddo
         r(i,i) = 1.d0
      enddo

c     # Only rotation about the z-axis, I believe.
      r(1,1) = cos(th)
      r(1,2) = sin(th)
      r(2,1) = -r(1,2)
      r(2,2) = r(1,1)

      do i = 1,3
         vrot(i) = wr(i)
      enddo

      do i = 1,3
         vrot(i) = 0
         do k = 1,3
            vrot(i) = vrot(i) + r(k,i)*wr(k)
         enddo
      enddo

      do i = 1,3
         wr(i) = vrot(i)
      enddo

      end
