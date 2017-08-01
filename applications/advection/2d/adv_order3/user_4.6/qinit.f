      subroutine clawpack46_qinit(maxmx,maxmy, meqn,mbc,mx,my,
     &      xlower,ylower,dx,dy,q,maux,aux)
       implicit none

       integer meqn, mbc, mx, my, maux, maxmx, maxmy
       double precision xlower, ylower, dx, dy, xi, yj
       double precision q(1-mbc:maxmx+mbc, 1-mbc:maxmy+mbc, meqn)
       double precision aux(1-mbc:maxmx+mbc, 1-mbc:maxmy+mbc, maux)
       double precision a(3,3), f, x, y

       integer i, j, mq, l, k
       

       f(x,y)=dsin(8d0*datan(1d0)*x)*dsin(8d0*datan(1d0)*y)


       a(1,1)=1d0
       a(1,2)=4d0
       a(1,3)=1d0
       a(2,1)=4d0
       a(2,2)=16d0
       a(2,3)=4d0
       a(3,1)=1d0
       a(3,2)=4d0
       a(3,3)=1d0

      do 20 j=1-mbc,my+mbc
         do 10 i=1-mbc,mx+mbc
            q(i,j,1)=0d0
   10       continue
   20    continue

      do mq=1, meqn 
         do j=1-mbc,my+mbc
            yj = ylower + (j-1d0)*dy
            do i=1-mbc,mx+mbc
               xi = xlower + (i-1d0)*dx
               do k=0,2
                  do l=0,2                  
                     q(i,j,mq)=q(i,j,mq)+(1d0/36d0)*a(k+1,l+1)
     &                     * f(dx*k/2d0+xi,dy*l/2d0+yj)
                  enddo
               enddo
            enddo
         enddo
       enddo

       return
       end
