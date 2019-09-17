subroutine tri_solve(N,dl,d,du,b)
    implicit none

    integer N
    double precision, dimension(N) :: dl, d, du, b

    integer i
    double precision fact

    do i = 1,N-2
        if (abs(d(i)) .ne. 0) then
            fact = dl(i)/d(i)
            d(i+1) = d(i+1) - fact*du(i)
            b(i+1) = b(i+1) - fact*b(i)       
        endif
        dl(i) = 0
    end do

111  format('b(',I3,') = ', 5E20.8)                 

    i = N-1
    if (abs(d(i)) .ne. 0) then
        fact = dl(i)/d(i)
        d(i+1) = d(i+1) - fact*du(i)
        b(i+1) = b(i+1) - fact*b(i)
    endif
    
    b(N) = b(N)/d(N)
    b(N-1) = (b(N-1) - du(N-1)*b(N))/d(N-1)
    do i = N-2,1,-1
        b(i) = (b(i) - du(i)*b(i+1) - dl(i)*b(i+2))/d(i)
    end do

end subroutine tri_solve