ccc   re-initialise incoming nucleon momenta in nucleon-nucleon cms
      subroutine pAinit(io)
      implicit none
      double precision sp,betaap,betap,betaa,beta
      integer io

      include 'mom.f'
      include 'vars.f'
      include 'varsi.f'
      include 'mp.f'
      include 'mion.f'
      include 'ion.f'
      include 'spA.f'
      include 'sAA.f'
      include 'ion_inel.f'

ccccc Used for incoherent photon emission from AA

      saa=an**2*s
      rtsaa=dsqrt(saa)

cccc   Using rigidity formula - ion energy per nucleon is Z/A x proton energy

      sp=si*(an/az)**2

      beta=dsqrt(1d0-4d0*mp**2/sp)
      betaa=dsqrt(1d0-4d0*mp**2/sp*an**2/az**2)

      spa=mp**2+mion**2+sp*az/2d0*(1d0+beta*betaa)
      rtspa=dsqrt(spa)
      betap=dsqrt(1d0-4d0*mp**2/spa)
      betaap=dsqrt(1d0-4d0*mion**2/spa)

      beta=dsqrt(1d0-4d0*mp**2/sp*an**2/az**2)
      betaa=dsqrt(1d0-4d0*mp**2/sp*an**2/az**2)
      spa=mp**2+mion**2+sp*az/2d0*(1d0+beta*betaa)*az/an
      rtspa=dsqrt(spa)
      
cccc  proton in +/- z direction  
      if(io.eq.1)then

      q(1,1)=0d0
      q(2,1)=0d0
      q(4,1)=rtspa/2d0-(mion**2-mp**2)/rtspa/2d0
      q(3,1)=dsqrt(q(4,1)**2-mp**2)

      q(1,2)=0d0
      q(2,2)=0d0
      q(4,2)=rtspa/2d0+(mion**2-mp**2)/rtspa/2d0
      q(3,2)=-dsqrt(q(4,2)**2-mion**2)

cccccccc     

      else

      q(1,1)=0d0
      q(2,1)=0d0
      q(4,1)=rtspa/2d0+(mion**2-mp**2)/rtspa/2d0
      q(3,1)=dsqrt(q(4,1)**2-mion**2)

      q(1,2)=0d0
      q(2,2)=0d0
      q(4,2)=rtspa/2d0-(mion**2-mp**2)/rtspa/2d0
      q(3,2)=-dsqrt(q(4,2)**2-mp**2)

      endif

      if(AA_frame)then

      if(io.eq.1)then

      q(1,1)=0d0
      q(2,1)=0d0
      q(4,1)=rtsnn/2d0
      q(3,1)=dsqrt(q(4,1)**2-mp**2)

      q(1,2)=0d0
      q(2,2)=0d0
      q(4,2)=an*rtsnn/2d0
      q(3,2)=-dsqrt(q(4,2)**2-mion**2)

      else

      q(1,1)=0d0
      q(2,1)=0d0
      q(4,1)=an*rtsnn/2d0
      q(3,1)=dsqrt(q(4,1)**2-mion**2)

      q(1,2)=0d0
      q(2,2)=0d0
      q(4,2)=rtsnn/2d0
      q(3,2)=-dsqrt(q(4,2)**2-mp**2)

      endif

      endif

      return
      end

