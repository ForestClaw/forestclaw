C     -*- FORTRAN -*-

c     # Headers files for C routines
      interface
         subroutine c_mult(s,n,x) bind(c,name="mult")
         use iso_c_binding, only : c_ptr
         implicit none

         type(c_ptr) s
         integer n
         double precision x(n)   !! double precision works too
         end subroutine
      end interface
