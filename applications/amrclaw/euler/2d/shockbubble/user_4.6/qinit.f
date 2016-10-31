      subroutine clawpack46_qinit(maxmx,maxmy,meqn,mbc,
     &      mx,my,xlower,ylower,dx,dy,q,maux,aux)
c     =====================================================
c
c     # Set initial conditions for q.
c
      implicit double precision (a-h,o-z)
      dimension q(1-mbc:maxmx+mbc, 1-mbc:maxmy+mbc, meqn)

       integer blockno, fc2d_clawpack46_get_block

       common /comic/ qin(5),qout(5)
       common /cominf/ rinf,vinf,einf
c
c
       blockno = fc2d_clawpack46_get_block()

       do 50 i=1-mbc,mx+mbc
          xclow = xlower + (i-1)*dx
          do 20 j=1-mbc,my+mbc
             yclow = ylower + (j-1)*dy
c            # set (xlow,ylow) to lower left corner of grid cell:
             call cellave2(blockno,xclow,yclow,dx,dy,win)
c            # win is now the fraction of the cell that lies inside the circle
             do 10 m=1,meqn
                q(i,j,m) = win*qin(m) + (1.d0-win)*qout(m)
   10        continue
   20     continue

          if (xclow .lt. 0.2d0) then
c            # behind shock:
             do 30 j=1-mbc,my+mbc
                q(i,j,1) = rinf
                q(i,j,2) = rinf*vinf
                q(i,j,3) = 0.d0
                q(i,j,4) = einf
                if (meqn .eq. 5) then
                   q(i,j,5) = 1.d0
                endif
   30        continue
          end if
c
          if (xclow .lt. 0.5d0 .and. meqn .eq. 5) then
c            # to give two different values of tracer in bubble
c            # for better visualization of motion:
             do 40 j=1,my
                q(i,j,5) = 2.d0*q(i,j,5)
   40        continue
          end if

   50  continue
       return
       end
