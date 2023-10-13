      function screen_01(qt)
      implicit none
      double precision screen_01,wt,sum,bt,hbt,qt,lbt,lbtmin,lbtmax,hlbt
      double precision opacpbint,opacpbint_3
      double precision b0,hbt0,sum0,btminb
      double precision btmin,btmax,aj0,aj0t,opac0,sumt
      integer iphi,ibt,nphi,nbt,nbt1,nbt2,nbt0
!if defined(__INTEL_COMPILER) || defined(__INTEL_LLVM_COMPILER)
      double precision dbesj0 
!endif
      include 'pi.f'
      include 'ion.f'
      include 'beam.f'
      include 'p0Xn.f'
      include 'opacpbpars.f'

      btmin=rzg*1.6d0

      nbt0=3000
      nbt=3000
      b0=4d0*rzg
      hbt0=(b0-btmin)/dble(nbt0)
      sum0=0d0
      
cccccccccccc

 111  btminb=0d0
      btmax=4d0*rzg

      aj0=(1d0-opacpbint_3(rzg*40d0))*rzg*40d0
      aj0=-aj0/2d0/pi

      hbt=(btmax-btminb)/dble(nbt)
      sumt=0d0

      do ibt=1,nbt

         bt=btminb+(dble(ibt)-0.5d0)*hbt

         wt=aj0

         opac0=1d0-opacpbint_3(bt)
         opac0=opac0/(-2d0*pi*aj0)*bt

         wt=wt*(opac0-1d0)
         wt=wt*dbesj0(qt*bt)
         wt=wt*hbt

         sumt=sumt+wt

      enddo

      sum=sumt

ccccc Finally add in analytic integration of J0 over 0->Infty

      screen_01=sum+aj0/qt
      
      return
      end
