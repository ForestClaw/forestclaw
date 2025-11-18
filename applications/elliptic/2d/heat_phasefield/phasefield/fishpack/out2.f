c
c
c =========================================================
      subroutine out2(maxmx,maxmy,mx,my,xlower,ylower,
     &                 dx,dy,phi,u,t,iframe)
c =========================================================
c
c     # Output the results for a general system of conservation laws
c     # in 2 dimensions
c
      implicit double precision (a-h,o-z)
      dimension phi(maxmx, maxmy)
      dimension u(maxmx, maxmy)
      character*10 fname1, fname2
      integer grid_number, amr_level,mpi_rank,block_number
      integer maux, ngrids
c
c     # Write the results to the file fort.q<iframe>
c     # Use format required by matlab script  plotclaw2.m
c     # The same format is used by the amrclaw package.
c     # Here it's adapted to output just the single grid.
c
c     # first create the file name and open file
c
      meqn = 2
      maux = 0
      ngrids = 1
      call header_ascii(iframe,t,meqn,maux,ngrids)


      fname1 = 'fort.qxxxx'
      nstp = iframe
      do 55 ipos = 10, 7, -1
         idigit = mod(nstp,10)
         fname1(ipos:ipos) = char(ichar('0') + idigit)
         nstp = nstp / 10
 55   continue

      open(unit=50,file=fname1,status='unknown',form='formatted')

c
c     # the following parameters are used in amrclaw where there are
c     # multiple grids.  Here they are all set to 1:

      grid_number = 0
      amr_level = 0
      block_number = 0
      mpi_rank = 0

      write(50,1001) grid_number,amr_level,block_number,mpi_rank,mx,my
 1001 format(I5,'                 grid_number', /,
     &       I5,'                 amr_level', /,
     &       I5,'                 block_number', /,
     &       I5,'                 mpi_rank',/,
     &       i5,'                 mx',/,
     &       i5,'                 my')

      write(50,1002) xlower,ylower,dx,dy
 1002 format(e18.8,'    xlow', /,
     &       e18.8,'    ylow', /,
     &       e18.8,'    dx', /,
     &       e18.8,'    dy',/)
c
      do 20 j=1,my
        do 10 i=1,mx
c          # exponents with more than 2 digits cause problems reading
c          # into matlab... reset tiny values to zero:
           if (dabs(phi(i,j)) .lt. 1d-99) phi(i,j) = 0.d0
           if (dabs(u(i,j)) .lt. 1d-99) u(i,j) = 0.d0
c
          write(50,1005) u(i,j), phi(i,j)
 1005     format(4e24.16)
c
 10       continue
        write(50,*) ' '
 20     continue
      write(50,*) ' '

      close(unit=50)

      return
      end


      subroutine header_ascii(iframe,time,meqn,maux,ngrids)
      implicit none

      integer iframe,meqn,ngrids, maux

      double precision time
      character*10 fname2
      integer matunit2, nstp, ipos, idigit

      matunit2 = 15

      fname2 = 'fort.txxxx'
      nstp = iframe
      do 55 ipos = 10, 7, -1
         idigit = mod(nstp,10)
         fname2(ipos:ipos) = char(ichar('0') + idigit)
         nstp = nstp / 10
 55   continue


      open(unit=matunit2,file=fname2)
      write(matunit2,1000) time,meqn,ngrids,maux,2
 1000 format(e30.20,'    time', /,
     &      i5,'                 meqn'/,
     &      i5,'                 ngrids'/,
     &      i5,'                 num_aux'/,
     &      i5,'                 num_dim')

      close(matunit2)

      end

