SUBROUTINE clawpack5_src2(meqn,mbc,mx,my,xlower,ylower, &
     dx,dy,q,maux,aux,t,dt)
  IMPLICIT NONE

  INTEGER meqn,mbc,mx,my,maux
  DOUBLE PRECISION xlower,ylower,dx,dy,t,dt
  DOUBLE PRECISION    q(meqn,1-mbc:mx+mbc, 1-mbc:my+mbc)
  DOUBLE PRECISION  aux(maux,1-mbc:mx+mbc, 1-mbc:my+mbc)

  DOUBLE PRECISION qstar(4)
  DOUBLE PRECISION gamma,gamma1
  COMMON /cparam/  gamma,gamma1

  INTEGER i,j,ndim
  DOUBLE PRECISION dt2,press,rad,rho,u,v

  !! # 2-stage Runge-Kutta method

  dt2    = dt/2.d0
  press  = 0.d0
  ndim   = 2

  DO  i = 1,mx
     DO j = 1,my
        rad      = aux(1,i,j)
        rho      = q(1,i,j)
        u        = q(2,i,j)/q(1,i,j)
        v        = q(3,i,j)/q(1,i,j)
        press    = gamma1*(q(4,i,j) - 0.5d0*rho*(u**2 + v**2))

        IF (rad.EQ.0.d0) WRITE(6,*) 'rad = 0 in src2'
        IF (rho.EQ.0.d0) WRITE(6,*) 'rho = 0 in q'

        qstar(1) = q(1,i,j) - dt2*(ndim-1)/rad * q(3,i,j)
        qstar(2) = q(2,i,j) - dt2*(ndim-1)/rad * (rho*u*v)
        qstar(3) = q(3,i,j) - dt2*(ndim-1)/rad * (rho*v*v)
        qstar(4) = q(4,i,j) - dt2*(ndim-1)/rad * v*(q(4,i,j) + press)

        !! # second stage

        rho      = qstar(1)
        u        = qstar(2)/qstar(1)
        v        = qstar(3)/qstar(1)
        press    = gamma1*(qstar(4) - 0.5d0*rho*(u**2 + v**2))
        IF (rho.EQ.0.d0) WRITE(6,*) 'rho = 0 in qstar'

        q(1,i,j) = q(1,i,j) - dt*(ndim-1)/rad * qstar(3)
        q(2,i,j) = q(2,i,j) - dt*(ndim-1)/rad * (rho*u*v)
        q(3,i,j) = q(3,i,j) - dt*(ndim-1)/rad * (rho*v*v)
        q(4,i,j) = q(4,i,j) - dt*(ndim-1)/rad * v*(qstar(4) + press)
     ENDDO
  ENDDO

  RETURN
END SUBROUTINE clawpack5_src2
