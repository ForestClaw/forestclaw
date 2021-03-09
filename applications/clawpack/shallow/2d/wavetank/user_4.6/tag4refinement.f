      subroutine wavetank_fort_tag4refinement(mx,my,mbc,
     &      meqn, xlower,ylower,dx,dy,time, blockno,
     &      q, tag_threshold, level, init_flag,tag_patch)
      implicit none

      integer mx,my, mbc, meqn, tag_patch, init_flag
      integer blockno, level
      double precision xlower, ylower, dx, dy, time
      double precision tag_threshold
      double precision q(1-mbc:mx+mbc,1-mbc:my+mbc,meqn)
      double precision xp,yp,zp

      integer*8 cont, get_context   
      logical fclaw2d_map_is_used

      integer i,j, mq, maxlevel
      double precision qmin, qmax, xc, yc, bathy, hij

      DOUBLE PRECISION xvals(5), bvals(5), xm, b


      DATA xvals /-160.8, -125.5, -90,    0, 40/
      DATA bvals /-4,     -4,      -0.45, 0,  2/

      tag_patch = 0

      cont = get_context()

      xm = -10

      maxlevel = 3



cc     # Do 1d refinement
c      do i = 1,mx
c         xc = xlower + (i-0.5)*dx
c         if (xc .gt. xvals(4)) then
c            return            
c         else if (xvals(3) .le. xc .and. xc .lt. xvals(4)) then
cc           # Refine to maxlevel            
c            tag_patch = 1
c            return
c         else if (xvals(2) .le. xc .and. xc .lt. xvals(3)) then
cc           # Refine to maxlevel - 1        
c            if (level .lt. maxlevel-1) then
c               tag_patch = 1
c               return
c            endif
c         endif
c      end do


c     !! Tag based on surface height      
      mq = 1
      qmin = 100
      qmax = -100
      do i = 1,mx
         xc = xlower + (i-0.5)*dx
         !! SEveral conditions where we shouldn't refine.
         if (time .ge. 15 .and. xc .lt.  -140) then
              return
         else if (time .ge. 30 .and. xc < -100) then            
            return
         elseif (time .ge. 45 .and. xc < -70) then
            return
         elseif (time .ge. 60 .and. xc < -50) then
             return
         endif
         b = bathy(xc,0)
         hij = q(i,j,1)
         qmin = min(hij + b,qmin)
         qmax = max(hij + b,qmax)
         if (hij .gt. 0 .and. abs(qmax) .gt. tag_threshold) then
            tag_patch = 1
            return
         endif
      enddo


      return



c     # Refine based only on first variable in system.
      mq = 1
      do j = 1,my
         do i = 1,mx
            xc = xlower + (i-0.5)*dx
            if (xm .le. xc .and. xc .le. xvals(4)+5) then
               tag_patch = 1
               return
            endif
         enddo
      enddo

      end
