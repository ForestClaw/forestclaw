      subroutine qinit_transport(mx,my,meqn,mbc,
     &      xlower,ylower,dx,dy,q,maux,aux,
     &      blockno, cont, xp,yp,zp)

      implicit none

      integer meqn, mbc, mx, my, maux,this_block_idx
      double precision xlower, ylower, dx, dy
      integer*8 cont
      integer blockno
      double precision q(1-mbc:mx+mbc, 1-mbc:my+mbc, meqn)
      double precision aux(1-mbc:mx+mbc, 1-mbc:my+mbc, maux)

      double precision xp(-mbc:mx+mbc+1,-mbc:my+mbc+1)
      double precision yp(-mbc:mx+mbc+1,-mbc:my+mbc+1)
      double precision zp(-mbc:mx+mbc+1,-mbc:my+mbc+1)

      integer i,j, ichoice, get_init_choice
      double precision x,y,z, xlow, ylow, w
      double precision gaussian_sum, cosine_bell_sum
      double precision correlated_bell,slotted_disk_sum

      ichoice = get_init_choice()

      call set_block(blockno)

      do j = 1-mbc,my+mbc
         do i = 1-mbc,mx+mbc
c            xlow = xlower + (i-1)*dx
c            ylow = ylower + (j-1)*dy
            x = xp(i,j)
            y = yp(i,j)
            z = zp(i,j)

            if (ichoice .eq. 1) then
               q(i,j,1) = gaussian_sum(x,y,z)
            elseif (ichoice .eq. 2 .or. ichoice .eq. 3) then
               q(i,j,1) = cosine_bell_sum(x,y,z)
               if (ichoice .eq. 3) then
                  q(i,j,2) = correlated_bell(q(i,j,1))
               endif
            elseif (ichoice .eq. 4) then
c               call cellave2(xlow,ylow,dx,dy,w)
               q(i,j,1) = slotted_disk_sum(x,y,z)
            endif
         enddo
      enddo

      return
      end


      double precision function  fdisc(xc,yc)
      implicit none

      double precision xc,yc, xp, yp, zp
      double precision q, slotted_disk_sum
      integer*8 cont, get_context
      integer blockno, get_block

      cont = get_context()
      blockno = get_block()

      call fclaw2d_map_c2m(cont,blockno,xc,yc,xp,yp,zp)

c     # call mapc2m(xc,yc,xp,yp,zp)

c     # Returns 0 or 1.
      q = slotted_disk_sum(xp,yp,zp)

      fdisc = q
      end



c     # --------------------------------------------------------
c     # This is to initialize q0 needed for exact solution.
c     # --------------------------------------------------------
      subroutine compute_q0(mx,my,mbc,xp,yp,zp,q0)
      implicit none

      integer mx,my,mbc
      double precision xp(-mbc:mx+mbc+1,-mbc:my+mbc+1)
      double precision yp(-mbc:mx+mbc+1,-mbc:my+mbc+1)
      double precision zp(-mbc:mx+mbc+1,-mbc:my+mbc+1)
      double precision q0(0:mx+1,0:my+1)

      integer i,j, ichoice, get_init_choice
      double precision x,y,z,t
      double precision gaussian_sum
      double precision cosine_bell_sum
      double precision slotted_disk_sum
      double precision w(3)

      ichoice = get_init_choice()

      t = 0
      do i = 0,mx+1
         do j = 0,my+1
            x = xp(i,j)
            y = yp(i,j)
            z = zp(i,j)
            if (ichoice .eq. 1) then
               q0(i,j) = gaussian_sum(x,y,z)
            elseif (ichoice .eq. 2 .or. ichoice .eq. 3) then
               q0(i,j) = cosine_bell_sum(x,y,z)
            elseif (ichoice .eq. 4) then
               q0(i,j) = slotted_disk_sum(x,y,z)
            endif
         enddo
      enddo

      end

      double precision function correlated_bell(q1)
      implicit none

      double precision q1

      double precision a,b
c     # Track the correlated cosine bell along with cosine bell
      a = -0.8d0
      b = 0.9d0

c     # Is this (b + cq)^2, from above?  Or just q^2
      correlated_bell = a*q1**2 + b

      return
      end
