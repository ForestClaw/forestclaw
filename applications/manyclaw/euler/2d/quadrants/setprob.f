c     ==================
      subroutine setprob
c     ==================

      implicit double precision (a-h,o-z)
      common /param/  gamma,gamma1
c
      open(unit=7,file='setprob.data',status='old',form='formatted')

c
       read(7,*) gamma
       gamma1 = gamma - 1.d0

       call set_maptype()

      return
      end
