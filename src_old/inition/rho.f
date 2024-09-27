      function rhoz(r)
      implicit none
      double precision rhoz,r

      include 'rho0.f'
      include 'ion.f'

      rhoz=rho0z/(1d0+dexp((r-Rzg)/dzg))

      return
      end

      function rhon(r)
      implicit none
      double precision rhon,r

      include 'rho0.f'
      include 'ion.f'

      rhon=rho0n/(1d0+dexp((r-Rng)/dng))

      return
      end
