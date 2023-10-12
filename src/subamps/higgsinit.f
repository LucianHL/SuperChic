ccc   gg --> SM higgs subprocess amplitude : normalisation
      subroutine higgsinit(mx)
      implicit none
      double precision xh,kf,ih,hmpp,fh,mx,alphas

      include 'pi.f'
      include 'ewpars.f'
      include 'normh.f'

      xh=(mt/mx)**2
      fh=-2d0*(dasin(1d0/(2d0*dsqrt(xh))))**2
      ih=3d0*xh*(2d0+(4d0*xh-1d0)*fh)

      hmpp=ih*mx**2*alphas(mx**2/4d0)/4d0/pi/v*2d0/3d0
      kf=(1d0+alphas(mx**2/4d0)/pi*(pi**2+11d0/2d0))

      normh=hmpp*dsqrt(kf)

      return
      end







