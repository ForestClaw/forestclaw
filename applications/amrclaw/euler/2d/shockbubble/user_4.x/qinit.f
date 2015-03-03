
c
c
c
c     =====================================================
       subroutine qinit(maxmx,maxmy,meqn,mbc,mx,my,xlower,ylower,
     &                   dx,dy,q,maux,aux)
c     =====================================================
c
c     # Set initial conditions for q.
c
       implicit double precision (a-h,o-z)
       dimension q(1-mbc:maxmx+mbc, 1-mbc:maxmy+mbc, meqn)

       integer blockno, fc2d_clawpack46_get_block

       integer*8 cont, get_context

       common /comic/ qin(5),qout(5)
       common /cominf/ rinf,vinf,einf
c
c
       blockno = fc2d_clawpack46_get_block()
       cont = get_context()

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

          call fclaw2d_map_c2m(cont,blockno,xclow,yclow,xlow,ylow,zlow)
          if (xlow .lt. 0.2d0) then
c            # behind shock:
             do 30 j=1-mbc,my+mbc
                q(i,j,1) = rinf
                q(i,j,2) = rinf*vinf
                q(i,j,3) = 0.d0
                q(i,j,4) = einf
                q(i,j,5) = 0.d0
   30        continue
          end if
c
c         if (xlow .lt. 0.5d0) then
c           # to give two different values of tracer in bubble
c           # for better visualization of motion:
c           do 40 j=1,my
c              q(i,j,5) = 2.d0*q(i,j,5)
c  40          continue
c           end if

   50    continue
       return
       end
