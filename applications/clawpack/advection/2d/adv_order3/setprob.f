      subroutine periodic_setprob()
      implicit none

      double precision pi,  pi2
      common /compi/ pi, pi2

      double precision ubar, vbar
      common /comrp/ ubar,vbar

      integer example
      common /example_comm/ example  

      integer mapping
      common /mapping_comm/ mapping

      pi = 4.d0*datan(1.d0)
      pi2 = 2.d0*pi


      open(unit=7,file='setprob.data',status='old',form='formatted')
      read(7,*) ubar
      read(7,*) vbar
      close(7)

      write(6,*) '-------------------> Setting example = mapping = 0'
      example = 0
      mapping = 0

      end
