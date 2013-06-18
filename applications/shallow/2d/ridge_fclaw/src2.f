c
c      =======================================================
       subroutine src2(maxmx,maxmy,meqn,mbc,mx,my,xlower,ylower,
     &                 dx,dy,q,maux,aux,t,dt)
c      =======================================================
c
       implicit double precision (a-h,o-z)
       dimension    q(1-mbc:maxmx+mbc, 1-mbc:maxmy+mbc, meqn)
       dimension aux(1-mbc:maxmx+mbc,1-mbc:maxmy+mbc, *)
       double precision RK(4,3)
       common /comsphere/ Rsphere, Omega
       common /sw/  g

c
c     # source term routine for Rossby-Haurwitz wave
c     # the source term models the Coriolis force using a 4-stage RK method
c     # and the projection of the velocity components to the tangent plane
c
c     df=12.600576e0 ! what is this??
      df = 2.d0*Omega

c     # project momentum components of q onto tangent plane:

      return

      do i=1,mx
        do j=1,my
            erx = aux(i,j,14)
            ery = aux(i,j,15)
            erz = aux(i,j,16)
            qn = erx*q(i,j,2) + ery*q(i,j,3) + erz*q(i,j,4)

            q(i,j,2) = q(i,j,2) - qn*erx
            q(i,j,3) = q(i,j,3) - qn*ery
            q(i,j,4) = q(i,j,4) - qn*erz

            enddo
        enddo

      return
c     ### no coriolis! ####################################

c     # calculate Coriolis term
      do i=1,mx
        xc = xlower + (i-0.5d0)*dx
        do j=1,my
            yc = ylower + (j-0.5d0)*dy
c
            call mapc2m(xc,yc,xp,yp,zp)
            erx = xp
            ery = yp
            erz = zp

c
            fcor = df*erz / (Rsphere**2.)

c           stage 1
            hu = q(i,j,2)
            hv = q(i,j,3)
            hw = q(i,j,4)

            RK(1,1) = dt*fcor*(erz*hv-ery*hw)
            RK(1,2) = dt*fcor*(erx*hw-erz*hu)
            RK(1,3) = dt*fcor*(ery*hu-erx*hv)

c           stage 2
            hu = q(i,j,2) + 0.5d0*RK(1,1)
            hv = q(i,j,3) + 0.5d0*RK(1,2)
            hw = q(i,j,4) + 0.5d0*RK(1,3)

            RK(2,1) = dt*fcor*(erz*hv-ery*hw)
            RK(2,2) = dt*fcor*(erx*hw-erz*hu)
            RK(2,3) = dt*fcor*(ery*hu-erx*hv)

c           stage 3
            hu = q(i,j,2) + 0.5d0*RK(2,1)
            hv = q(i,j,3) + 0.5d0*RK(2,2)
            hw = q(i,j,4) + 0.5d0*RK(2,3)

            RK(3,1) = dt*fcor*(erz*hv-ery*hw)
            RK(3,2) = dt*fcor*(erx*hw-erz*hu)
            RK(3,3) = dt*fcor*(ery*hu-erx*hv)

c           stage 4
            hu = q(i,j,2) + 0.5d0*RK(3,1)
            hv = q(i,j,3) + 0.5d0*RK(3,2)
            hw = q(i,j,4) + 0.5d0*RK(3,3)

            RK(4,1) = dt*fcor*(erz*hv-ery*hw)
            RK(4,2) = dt*fcor*(erx*hw-erz*hu)
            RK(4,3) = dt*fcor*(ery*hu-erx*hv)

            do m=2,meqn
               q(i,j,m) = q(i,j,m)
     &                 + (RK(1,m-1) + 2.d0*RK(2,m-1)+
     &                 2.d0*RK(3,m-1) + RK(4,m-1))/6.d0
            enddo

        enddo
      enddo

c     # project momentum components of q onto tangent plane:

      do i=1-mbc,mx+mbc
        do j=1-mbc, my+mbc
            erx = aux(i,j,14)
            ery = aux(i,j,15)
            erz = aux(i,j,16)
            qn = erx*q(i,j,2) + ery*q(i,j,3) + erz*q(i,j,4)
c
            q(i,j,2) = q(i,j,2) - qn*erx
            q(i,j,3) = q(i,j,3) - qn*ery
            q(i,j,4) = q(i,j,4) - qn*erz

            enddo
        enddo


       return
       end
