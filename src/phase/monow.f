ccc   generates rho meson invariant mass distribution according to modified BW 
      subroutine bwmono(mout,jrho)
      implicit none
      double precision mout,jrho

      include 'mres.f'
      include 'pi.f'
      include 'monopar.f'
      
      jrho=mres*gamm/((mres**2-mout**2)**2+mres**2*
     &     gamm**2)/pi

      return
      end
