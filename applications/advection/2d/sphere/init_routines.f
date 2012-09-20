      subroutine set_init_choice(ichoice)
      implicit none

      integer ichoice, ichoice_com

      common /cominit/ ichoice_com

      ichoice_com = ichoice


      end

      integer function get_init_choice()
      implicit none

      common /cominit/ ichoice_com
      integer ichoice_com

      get_init_choice = ichoice_com

      end


      subroutine set_locations(w,n)
      implicit none

      integer i, n, m
      double precision w(3,n)

      double precision wc_com(3,100)
      integer n_com

      common /comexp1/ wc_com
      common /comexp2/ n_com

c     # Locations of initial exponentials
      do i = 1,n
         do m = 1,3
            wc_com(m,i) = w(m,i)
         enddo
      enddo
      n_com = n
      end

      subroutine get_locations(i,w)
      implicit none

      integer i,m
      double precision w(3)

      double precision wc_com(3,100)
      integer n_com
      common /comexp1/ wc_com
      common /comexp2/ n_com

      if (i .gt. n_com) then
         write(6,*) 'get_locations : i .gt. n_com'
         stop
      endif

      do m = 1,3
         w(m) = wc_com(m,i)
      enddo

      end

      integer function n_locations()
      implicit none

      integer n_com
      common /comexp2/ n_com

      n_locations = n_com


      end


      subroutine set_wind_parms(kappa,tfinal)
      implicit none

      double precision kappa,tfinal
      double precision init_kappa_com, tfinal_com
      common /comwind/ init_kappa_com, tfinal_com

      init_kappa_com = kappa
      tfinal_com = tfinal

      end


      subroutine get_wind_parms(kappa,tfinal)
      implicit none

      double precision kappa,tfinal
      double precision init_kappa_com, tfinal_com
      common /comwind/ init_kappa_com, tfinal_com


      kappa = init_kappa_com
      tfinal = tfinal_com

      end
