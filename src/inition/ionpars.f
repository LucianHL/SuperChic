      subroutine ionpars
      implicit none
      integer fit
      
      include 'ion.f'
      include 'pi.f'
      
      rzg=rz/0.1973d0
      dzg=dz/0.1973d0
      rng=rn/0.1973d0
      dng=dn/0.1973d0

      rcore=0.8d0/0.1973d0
      nshell=20d0
      
      return
      end


      subroutine opacpars(rts)
      implicit none
      double precision rts,rtst
      
      include 'onechannel.f'
   
      rtst=rts/1d3

ccc   Fit according to one-channel calculation
      
      sig0=122.9d0+4.88d0*rtst-0.056d0*rtst**2
      a=0.158d0+0.005d0*rtst-0.000138d0*rtst**2
      b=9.78d0+2.47d0*rtst-0.054d0*rtst**2
      c=0.438d0-0.00524d0*rtst+0.000108d0*rtst**2 
      sig0=sig0/0.389d0

      return
      end


      
