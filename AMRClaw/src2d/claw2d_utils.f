      subroutine inputparms(initial_dt, tfinal, max_cfl, nout, src_term,
     &      mcapa, maux,meqn,mwaves,maxmwaves,mthlim,mbc,mthbc,order)
      implicit none

      double precision initial_dt, tfinal, max_cfl
      integer nout, src_term, mcapa, maux, meqn, mwaves
      integer maxmwaves,mbc
      integer mthbc(4),mthlim(maxmwaves), order(2)

      integer mw, m


      open(55,file='claw2ez.data')

c     timestepping variables
      read(55,*) initial_dt
      read(55,*) tfinal
      read(55,*) max_cfl
      read(55,*) nout

      read(55,*) src_term
      read(55,*) mcapa
      read(55,*) maux

      read(55,*) meqn
      read(55,*) mwaves

      if (mwaves .gt. maxmwaves) then
         write(6,*) 'ERROR : (claw_utils.f) mwaves > maxmwaves'
         stop
      endif

      read(55,*) (mthlim(mw), mw=1,mwaves)

      read(55,*) mbc
      read(55,*) mthbc(1)
      read(55,*) mthbc(2)
      read(55,*) mthbc(3)
      read(55,*) mthbc(4)

      read(55,*) (order(m),m=1,2)

      close(55)
      end
