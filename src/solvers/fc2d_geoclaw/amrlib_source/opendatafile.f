
      subroutine opendatafile(iunit, fname)
c     #
!> Open the file fname and determine how many leading lines are
!! comments.  Then rewind and skip over the comment lines so that the
!! file is ready for reading data from.
!!
!! All comment lines must start with # in the first column.

      integer iunit, line, commentlines


c     This allows for a variable name length and obsolves (I think) the
c     length problem that appeared for some compilers (notably ifort)
      character(len=*) :: fname
c     character*25 fname
c     character*12 fname12

      character*1 firstchar
      logical foundFile
      
c     # Debug statement to test fname:
c     write(6,*) '+++ fname = XXX',fname,'XXX'

      inquire(file=fname,exist=foundFile)
      if (.not. foundFile) then
c         # for backward compatability: maybe fname is declared character*12
c         # in calling routine...
c         fname12 = fname(1:12)
c         #write(6,*) 'truncated fname12 = XXX',fname12,'XXX'
c         inquire(file=fname12,exist=foundFile)
          if (.not. foundFile) then
            print "(2a)",'*** in opendatafile, file not found:', fname 
            stop
          endif
c         open(unit=iunit,file=fname12,status='old',form='formatted')
c         write(6,*) 'Reading data file: ', fname12
      else
          open(unit=iunit,file=fname,status='old',form='formatted')
          print "(2a)",'Reading data file: ', fname
      endif

c     # this version may not work in f77
c     firstchar = '#'
c     commentlines = -1
c     do while (firstchar == '#')
c        read(iunit,*) firstchar
c        commentlines = commentlines + 1
c        enddo

      firstchar = '#'
      commentlines = -1
      do commentlines=-1,1000
          if (firstchar .eq. '#') then
              read(iunit,"(a1)") firstchar
          else
              exit
          endif
      enddo
      
      write(6,602) commentlines
  602 format('         first',i2,
     &       ' lines are comments and will be skipped')
     
c     Rewind file and go back to data line
C       print *,firstchar
C       rewind(iunit)
C       do line=1,commentlines
C          read(iunit,*) firstchar
C          print *,firstchar
C          enddo
C       
      return
      end
