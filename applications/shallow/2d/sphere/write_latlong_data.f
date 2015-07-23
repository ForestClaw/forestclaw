      subroutine write_latlong_data(lat, longitude)
      implicit none

      double precision lat(2), longitude(2)

c     # This file will be read by Matlab.
      open(10,file='latlong.dat')
      write(10,100) lat(1),lat(2),longitude(1),
     &      longitude(2)
      close(10)

  100 format(4E24.16)


      end
