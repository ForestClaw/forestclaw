SUBROUTINE clawpack46_setaux(maxmx,maxmy,mbc,mx,my, & 
           xlower,ylower,dx,dy,maux,aux)

    IMPLICIT NONE
    INTEGER maxmx, maxmy, mbc,mx,my, maux
    DOUBLE PRECISION xlower,ylower, dx, dy
    DOUBLE PRECISION  aux(1-mbc:mx+mbc,1-mbc:my+mbc,maux)

    INTEGER i, j, ibc, jbc
    DOUBLE PRECISION xc, bathy_compete, b, slope, d2xzb

    write(6,*) 'Calling setaux ...'

    DO i = 1-mbc,mx+mbc
        DO j = 1-mbc,my+mbc
            xc = xlower + (i-0.5)*dx     
            call bathy_complete(xc,b,slope,d2xzb)
            aux(i,j,1) = b
            !! aux(2,i) = slope
            !! aux(3,i) = d2xzb
        end do
    end do

    do i = 2-mbc,mx+1
        aux(i,2) = (aux(i+1,1) - aux(i-1,1))/(2.d0*dx)
        aux(i,3) = (aux(1,i+1) - 2.d0*aux(1,i) + aux(i-1,1))/(dx*dx)
    END DO



END SUBROUTINE CLAWPACK46_SETAUX