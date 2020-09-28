      subroutine periodic_setprob()
      implicit none

      double precision pi,  pi2
      common /compi/ pi, pi2

      integer example
      common /example_comm/ example  

      integer refinement_strategy
      common /refinement_comm/ refinement_strategy

      double precision ubar, vbar
      common /comrp/ ubar,vbar

      pi = 4.d0*datan(1.d0)
      pi2 = 2.d0*pi


      open(unit=7,file='setprob.data',status='old',form='formatted')
      read(7,*) example
      read(7,*) refinement_strategy
      read(7,*) ubar
      read(7,*) vbar
      close(7)
      end
