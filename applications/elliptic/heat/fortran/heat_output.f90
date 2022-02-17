subroutine heat_fort_header_ascii(matname1,matname2, time,meqn,maux,ngrids)
    implicit none

    integer meqn,ngrids, maux

    character*11 matname1
    character*11 matname2
    double precision time
    integer matunit1, matunit2
    integer mfields

    matunit1 = 10
    matunit2 = 15

    open(unit=matunit2,file=matname2)

    mfields = meqn + 2  !! include soln, error
    write(matunit2,1000) time,mfields,ngrids,maux,2
 1000 format(e30.20,'    time', /,         & 
           i5,'                 mfields'/,    & 
           i5,'                 ngrids'/,  & 
           i5,'                 num_aux'/, & 
           i5,'                 num_dim')

    close(matunit2)

    open(unit=matunit1,file=matname1,status='replace')
    close(matunit1)

end subroutine heat_fort_header_ascii

subroutine heat_fort_output_ascii(matname1, & 
         mx,my,meqn,mbc, xlower,ylower, dx,dy,  & 
         q,soln,error,patch_num,level,blockno,mpirank)

    implicit none

    character(len=11) matname1
    integer meqn,mbc,mx,my
    integer patch_num
    integer level, blockno, mpirank
    double precision xlower, ylower,dx,dy

    double precision q(1-mbc:mx+mbc,1-mbc:my+mbc,meqn)
    double precision error(1-mbc:mx+mbc,1-mbc:my+mbc,meqn)
    double precision soln(1-mbc:mx+mbc,1-mbc:my+mbc,meqn)

    double precision qvec(3), dmax, heat_eval_refinement
    integer matunit1
    integer i,j,mq

    matunit1 = 10
    open(matunit1,file=matname1,position='append');

    call fclaw2d_clawpatch46_fort_write_grid_header(matunit1, & 
           mx,my,xlower,ylower, dx,dy,patch_num,level, & 
           blockno,mpirank)

    if (meqn .gt. 5) then
        write(6,'(A,A,A,I5,A)')     & 
             'Warning (fclaw2d_fort_write_grid_header.f) ',  & 
             ': meqn > 5; change format statement 120.',     & 
             '(meqn = ',meqn,')'
        stop
    endif

    do j = 1,my
        do i = 1,mx
            do mq = 1,meqn
                if (abs(q(i,j,mq)) .lt. 1d-99) then
                    q(i,j,mq) = 0.d0
                elseif (abs(q(i,j,mq)) .gt. 1d99) then
                    q(i,j,mq) = 1d99
                endif
            end do
!!            if (abs(error(i,j,1)) .lt. 1d-99) then
!!                error(i,j,1) = 0.d0
!!            elseif (abs(error(i,j,1)) .gt. 1d99) then
!!                error(i,j,1) = 1d99
!!            if (abs(soln(i,j,1)) .lt. 1d-99) then
!!                soln(i,j,1) = 0.d0
!!            elseif (abs(soln(i,j,1)) .gt. 1d99) then
!!                soln(i,j,1) = 1d99
!!            endif
!!            write(matunit1,120) (q(i,j,mq),mq=1,meqn), qlap, error(i,j,1)
!!            endif

            qvec(1) = q(i,j,1)
            qvec(2) = q(i-1,j,1)
            qvec(3) = q(i,j-1,1)

            dmax = heat_eval_refinement(qvec,dx,dy)

            write(matunit1,120) (q(i,j,mq),mq=1,meqn), dmax
        end do
        write(matunit1,*) ' '
    end do

    close(matunit1)

120 format (5E26.16)
!!121 format (3I5,6E24.16)
!!122 format (2I5,6E24.16)

end subroutine heat_fort_output_ascii
