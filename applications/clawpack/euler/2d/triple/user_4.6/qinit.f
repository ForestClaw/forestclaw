      subroutine clawpack46_qinit(maxmx,maxmy,meqn,mbc,
     &      mx,my,xlower,ylower,dx,dy,q,maux,aux)
c     =====================================================
c
c     # Set initial conditions for q.
c
      implicit double precision (a-h,o-z)
      dimension q(1-mbc:maxmx+mbc, 1-mbc:maxmy+mbc, meqn)

      integer i, j, blockno, fc2d_clawpack46_get_block
      double precision xi, yj, win, qin(5), qout(5)
      double precision r1_xlower, r1_xupper
      double precision r1_ylower, r1_yupper
      double precision r2_xlower, r2_xupper
      double precision r2_ylower, r2_yupper
      double precision r3_xlower, r3_xupper
      double precision r3_ylower, r3_yupper

      common /comic/ qin, qout
      common /cparam/  gamma, gamma1

      blockno = fc2d_clawpack46_get_block()

      ! Region 1 Boundaries
      r1_xlower = 0
      r1_xupper = 1
      r1_ylower = -1.5
      r1_yupper = 1.5

      ! Region 2 Boundaries
      r2_xlower = 1
      r2_xupper = 7
      r2_ylower = 0
      r2_yupper = 1.5

      ! Region 3 Boundaries
      r3_xlower = 1
      r3_xupper = 7
      r3_ylower = -1.5
      r3_yupper = 0

      do i = 1-mbc,mx+mbc
         xi = xlower + (i - 0.5d0)*dx
         do j = 1-mbc,my+mbc
            yj = ylower + (j - 0.5d0)*dy

            call cellave2(blockno, xi, yj, dx, dy, win)
            do m = 1,meqn
               q(i,j,m) = win*qin(m) + (1.0d0 - win)*qout(m)
            enddo

            ! Region 1
            if (xi .gt. r1_xlower .and. xi .lt. r1_xupper) then
               if (yj .gt. r1_ylower .and. yj .lt. r1_yupper) then
                  q(i,j,1) = 1.0d0 ! Density
                  q(i,j,2) = 0.0d0 ! x-momentum
                  q(i,j,3) = 0.0d0 ! y-momentum
                  q(i,j,4) = 1.0d0 / gamma ! Energy
                  if (meqn .eq. 5) then
                     q(i,j,5) = 1.0d0 ! Tracer
                  endif
               endif
            endif

            ! Region 2
            if (xi .gt. r2_xlower .and. xi .lt. r2_xupper) then
               if (yj .gt. r2_ylower .and. yj .lt. r2_yupper) then
                  q(i,j,1) = 0.125d0 ! Density
                  q(i,j,2) = 0.0d0 ! x-momentum
                  q(i,j,3) = 0.0d0 ! y-momentum
                  q(i,j,4) = 0.1d0 / gamma ! Energy
                  if (meqn .eq. 5) then
                     q(i,j,5) = 2.0d0 ! Tracer
                  endif
               endif
            endif

            ! Region 3
            if (xi .gt. r3_xlower .and. xi .lt. r3_xupper) then
               if (yj .gt. r3_ylower .and. yj .lt. r3_yupper) then
                  q(i,j,1) = 1.0d0 ! Density
                  q(i,j,2) = 0.0d0 ! x-momentum
                  q(i,j,3) = 0.0d0 ! y-momentum
                  q(i,j,4) = 0.1d0 / gamma ! Energy
                  if (meqn .eq. 5) then
                     q(i,j,5) = 3.0d0 ! Tracer
                  endif
               endif
            endif

         enddo
      enddo

      !  idisc_example = idisc
!        do  i = 1-mbc,mx+mbc
!           xclow = xlower + (i-1)*dx
!           do j = 1-mbc,my+mbc
!              yclow = ylower + (j-1)*dy
! c            # set (xlow,ylow) to lower left corner of grid cell:

!              if (idisc_example .eq. 3) then
!                 idisc = 3
!                 call cellave2(blockno,xclow,yclow,dx,dy,win)
! c               # win is now the fraction of the cell that lies inside the
! c               # circle
!                 do m = 1,meqn
!                    q(i,j,m) = win*qin(m) + (1.d0-win)*qout(m)
!                 enddo
!                 idisc = 4
!                 call cellave2(blockno,xclow,yclow,dx,dy,win)
!                 q(i,j,5) = q(i,j,5) + win
!              else
!                 call cellave2(blockno,xclow,yclow,dx,dy,win)
! c               # win is now the fraction of the cell that lies inside the
! c               # circle
!                 do m=1,meqn
!                    q(i,j,m) = win*qin(m) + (1.d0-win)*qout(m)
!                 enddo
!              endif
!           enddo

! c         # behind shock:
!           do  j=1-mbc,my+mbc
!              if (xclow .lt. 0.2d0) then
!                 q(i,j,1) = rinf
!                 q(i,j,2) = rinf*vinf
!                 q(i,j,3) = 0.d0
!                 q(i,j,4) = einf
!                 if (meqn .eq. 5) then
!                    q(i,j,5) = 1.d0
!                 endif
!              endif
!           enddo
! c
!           if (xclow .lt. 0.5d0 .and. meqn .eq. 5) then
! c            # to give two different values of tracer in bubble
! c            # for better visualization of motion:
!              do j=1,my
!                 q(i,j,5) = 2.d0*q(i,j,5)
!              enddo
         !  end if
!        enddo

       return
       end
