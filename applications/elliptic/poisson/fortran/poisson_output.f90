subroutine poisson_fort_header_ascii(matname1,matname2, time,mfields,maux,ngrids)
    implicit none

    integer mfields,ngrids, maux

    character*11 matname1
    character*11 matname2
    double precision time
    integer matunit1, matunit2
    integer mfieldsp2

    matunit1 = 10
    matunit2 = 15

    open(unit=matunit2,file=matname2)

    mfieldsp2 = mfields + 2  !! include soln, error
    write(matunit2,1000) time,mfieldsp2,ngrids,maux,2
 1000 format(e30.20,'    time', /,         & 
           i5,'                 mfields'/,    & 
           i5,'                 ngrids'/,  & 
           i5,'                 num_aux'/, & 
           i5,'                 num_dim')

    close(matunit2)

    open(unit=matunit1,file=matname1,status='replace')
    close(matunit1)

end subroutine poisson_fort_header_ascii

subroutine poisson_fort_output_ascii(matname1, & 
         mx,my,mfields,mbc, xlower,ylower, dx,dy,  & 
         rhs,soln,error,patch_num,level,blockno,mpirank)

    implicit none

    character(len=11) matname1
    integer mfields,mbc,mx,my
    integer patch_num
    integer level, blockno, mpirank
    double precision xlower, ylower,dx,dy

    double precision rhs(1-mbc:mx+mbc,1-mbc:my+mbc,mfields)
    double precision error(1-mbc:mx+mbc,1-mbc:my+mbc,mfields)
    double precision soln(1-mbc:mx+mbc,1-mbc:my+mbc,mfields)

    integer matunit1
    integer i,j,mq

    matunit1 = 10
    open(matunit1,file=matname1,position='append');

    call fclaw2d_clawpatch46_fort_write_grid_header(matunit1, & 
           mx,my,xlower,ylower, dx,dy,patch_num,level, & 
           blockno,mpirank)

    if (mfields .gt. 5) then
        write(6,'(A,A,A,I5,A)')     & 
             'Warning (fclaw2d_fort_write_grid_header.f) ',  & 
             ': meqn > 5; change format statement 120.',     & 
             '(meqn = ',mfields,')'
        stop
    endif

    do j = 1,my
        do i = 1,mx
            do mq = 1,mfields
                if (abs(rhs(i,j,mq)) .lt. 1d-99) then
                    rhs(i,j,mq) = 0.d0
                elseif (abs(rhs(i,j,mq)) .gt. 1d99) then
                    rhs(i,j,mq) = 1d99
                endif
            end do
            if (abs(error(i,j,1)) .lt. 1d-99) then
                error(i,j,1) = 0.d0
            elseif (abs(error(i,j,1)) .gt. 1d99) then
                error(i,j,1) = 1d99
            endif
            if (abs(soln(i,j,1)) .lt. 1d-99) then
                soln(i,j,1) = 0.d0
            elseif (abs(soln(i,j,1)) .gt. 1d99) then
                soln(i,j,1) = 1d99
            endif
            write(matunit1,120) (rhs(i,j,mq),mq=1,mfields), soln(i,j,1), error(i,j,1)
        end do
        write(matunit1,*) ' '
    end do
120 format (5E26.16)

    close(matunit1)

end subroutine poisson_fort_output_ascii
