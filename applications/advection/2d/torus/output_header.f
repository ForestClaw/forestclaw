      subroutine torus_fort_write_header(matname1,matname2,
     &      time,meqn,ngrids)
      implicit none

      integer iframe,meqn,ngrids

      character*10 matname1
      character*10 matname2
      double precision time
      integer matunit1, matunit2
      integer mfields

c     # Just to make sure this file gets replaced,
c     # and not appended.
      matunit1 = 10
      open(unit=matunit1,file=matname1,status='replace')
      close(matunit1)

      matunit2 = 10
      open(unit=matunit2,file=matname2)

c     # Write out error as an extra field (why not just do this in C?)
      mfields = meqn + 2
      write(matunit2,1000) time,mfields,ngrids
 1000 format(e30.20,'    time', /,
     &      i5,'                 mfields'/,
     &      i5,'                 ngrids')

      close(matunit2)

      end
