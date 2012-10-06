      subroutine qinit_mapped(mx,my,meqn,mbc,xlower, ylower, dx, dy,
     &      xp, yp, zp, q_claw,maux,aux,blockno)
      implicit none

      integer mx,my,meqn,mbc, maux, blockno
      double precision xlower, ylower, dx, dy
      double precision q_claw(1-mbc:mx+mbc,1-mbc:my+mbc,meqn)
      double precision    aux(1-mbc:mx+mbc,1-mbc:my+mbc,maux)

      double precision xp(-mbc:mx+mbc+1,-mbc:my+mbc+1)
      double precision yp(-mbc:mx+mbc+1,-mbc:my+mbc+1)
      double precision zp(-mbc:mx+mbc+1,-mbc:my+mbc+1)

      double precision a,b,c, q0, q1,q2
      integer i,j, ichoice, get_init_choice
      double precision t,x,y,z
      double precision gaussian_sum, cosine_bell_sum
      double precision slotted_disk_sum

      ichoice = get_init_choice()

      if (ichoice .eq. 3 .and. meqn .eq. 1) then
         write(6,*) 'Error : Set meqn=2 for correlated cosine bells'
         stop
      endif

      a = -0.8d0
      b = 0.9d0

      t = 0
      do j = 1-mbc,my+mbc
         do i = 1-mbc,mx+mbc
            x = xp(i,j)
            y = yp(i,j)
            z = zp(i,j)

            if (x+y .le. 0) then
               q_claw(i,j,1) = 1.d0
            else
               q_claw(i,j,1) = 0.d0
            endif

c            if (ichoice .eq. 1) then
c               q1 = gaussian_sum(x,y,z)
c            elseif (ichoice .eq. 2 .or. ichoice .eq. 3) then
c               q1 = cosine_bell_sum(x,y,z)
c            elseif (ichoice .eq. 4) then
c               q1 = slotted_disk_sum(x,y,z)
c            endif
c            q_claw(i,j,1) = q1

            if (ichoice .eq. 3) then
               q2 = a*q1**2 + b
               q_claw(i,j,2) = q2
            endif
         enddo
      enddo


      end
