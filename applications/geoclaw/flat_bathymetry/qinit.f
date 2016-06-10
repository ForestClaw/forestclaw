! c     =====================================================
!        subroutine preqinit(min_level, max_level)
! c     =====================================================
!        use qinit_module

!        integer :: min_level
!        integer :: max_level

!        x_low_qinit = -9.0
!        y_low_qinit = -100.0
!        x_hi_qinit = 9.0
!        y_hi_qinit = 100.0

!        min_level_qinit = min_level
!        max_level_qinit = max_level
!        end
! c
! c
c
c     =====================================================
       subroutine qinit(meqn,mbc,mx,my,xlower,ylower,
     &                   dx,dy,q,maux,aux)
c     =====================================================
c
c     # Set initial conditions for q.
c     # Sample scalar equation with data that is piecewise constant with
c     # q = 1.0  if  0.1 < x < 0.6   and   0.1 < y < 0.6
c     #     0.1  otherwise
c
       use geoclaw_module, only: sea_level, coordinate_system
       use geoclaw_module, only: grav


       implicit none
       integer i,j,meqn,mbc,mx,my,maux
       double precision xlower,ylower,dx,dy,xi,yj
       double precision q(meqn,1-mbc:mx+mbc, 1-mbc:my+mbc)
       double precision aux(maux,1-mbc:mx+mbc, 1-mbc:my+mbc)
c
       q = 0.d0
       do 20 i=1-mbc,mx+mbc
c         # (xi, yj) is the locaiton of the physical domain
          xi = xlower + (i-0.5d0)*dx
          do 20 j=1-mbc,my+mbc
             yj = ylower + (j-0.5d0)*dy
             if (xi.lt.9.0 .and. xi.gt. -9.0) then
                q(1,i,j) = sea_level - aux(1,i,j) + 
     &           exp(-xi**2/18.0)
c     !            q(2,i,j) = sqrt(grav*100)* (-aux(1,i,j)+ 
c     ! &                        exp(-xi**2/18.0))
             else
                q(1,i,j) = sea_level - aux(1,i,j)
             endif
  20         continue
       return
       end
