      subroutine clawpack46_qinit(maxmx,maxmy,meqn,mbc,
     &      mx,my,xlower,ylower,dx,dy,q,maux,aux)
c     =====================================================
c
c     # Set initial conditions for q.
c
      implicit double precision (a-h,o-z)
      dimension q(1-mbc:maxmx+mbc, 1-mbc:maxmy+mbc, meqn)

       integer blockno, fc2d_clawpack46_get_block
       integer i,j, m
       double precision xclow, yclow, qin(5), qout(5)
       double precision rinf,vinf,einf, win, r0,x0,y0,alf,beta
       integer idisc, idisc_example

       common /comic/ qin,qout
       common /cominf/ rinf,vinf,einf
       common/cdisc/ x0,y0,alf,beta,r0,idisc
c
c
       blockno = fc2d_clawpack46_get_block()

       idisc_example = idisc
       do  i = 1-mbc,mx+mbc
          xclow = xlower + (i-1)*dx
          do j = 1-mbc,my+mbc
             yclow = ylower + (j-1)*dy
c            # set (xlow,ylow) to lower left corner of grid cell:

             if (idisc_example .eq. 3) then
                idisc = 3
                call cellave2(blockno,xclow,yclow,dx,dy,win)
c               # win is now the fraction of the cell that lies inside the
c               # circle
                do m = 1,meqn
                   q(i,j,m) = win*qin(m) + (1.d0-win)*qout(m)
                enddo
                idisc = 4
                call cellave2(blockno,xclow,yclow,dx,dy,win)
                q(i,j,5) = q(i,j,5) + win
             else
                call cellave2(blockno,xclow,yclow,dx,dy,win)
c               # win is now the fraction of the cell that lies inside the
c               # circle
                do m=1,meqn
                   q(i,j,m) = win*qin(m) + (1.d0-win)*qout(m)
                enddo
             endif
          enddo

c         # behind shock:
          do  j=1-mbc,my+mbc
             if (xclow .lt. 0.2d0) then
                q(i,j,1) = rinf
                q(i,j,2) = rinf*vinf
                q(i,j,3) = 0.d0
                q(i,j,4) = einf
                if (meqn .eq. 5) then
                   q(i,j,5) = 1.d0
                endif
             endif
          enddo
c
          if (xclow .lt. 0.5d0 .and. meqn .eq. 5) then
c            # to give two different values of tracer in bubble
c            # for better visualization of motion:
             do j=1,my
                q(i,j,5) = 2.d0*q(i,j,5)
             enddo
          end if
       enddo

       return
       end
