      function screen(qt)
      implicit none
      double precision screen,wt,sum,bt,hbt,qt
      double precision opacpbpint,opacpbint
      double precision b0,hbt0,sum0,btminb
      integer ibt,nbt,nbt0
!if defined(__INTEL_COMPILER) || defined(__INTEL_LLVM_COMPILER)
      double precision dbesj0 
!endif
      include 'pi.f'
      include 'ion.f'
      include 'beam.f'
      include 'p0Xn.f'
      include 'btmaxmin.f'
      include 'opacpbpars.f'
      include 'opacpbppars.f'

      if(beam.eq.'ion')nbt=ioppb
      if(beam.eq.'ionp')nbt=ioppbp

      nbt=3000

      if(ionbreakup)then
         if(fAA.eq.'00')then
            btmin=0d0            
         endif
         if(fAA.eq.'01')then
            nbt=3000
         endif
      else
         btmin=0d0
      endif

      goto 888

      nbt0=500
      nbt=8000
      b0=4d0*rzg
      hbt0=(b0-btmin)/dble(nbt0)
      sum0=0d0
      
      do ibt=1,nbt0

         bt=btmin+(dble(ibt)-0.5d0)*hbt0

         if(beam.eq.'ion')wt=1d0-opacpbint(bt)
         if(beam.eq.'ionp')wt=1d0-opacpbpint(bt)

         wt=wt*dbesj0(qt*bt)
         wt=-wt*bt*hbt0/2d0/pi

         sum0=sum0+wt

      enddo

 888  btminb=btmin

      hbt=(btmax-btminb)/dble(nbt)
      sum=0d0
      
      do ibt=1,nbt

         bt=btminb+(dble(ibt)-0.5d0)*hbt

         if(beam.eq.'ion')wt=1d0-opacpbint(bt)
         if(beam.eq.'ionp')wt=1d0-opacpbpint(bt)
 
         wt=wt*dbesj0(qt*bt)
         wt=-wt*bt*hbt/2d0/pi

         sum=sum+wt
         
      enddo

      sum=sum+sum0

      screen=sum
      
      return
      end
