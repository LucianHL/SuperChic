      function screen(qt)
      implicit none
      double precision screen,wt,sum,bt,hbt,qt,lbt,lbtmin,lbtmax,hlbt
      double precision opacpbpint,opacpbint
      double precision b0,hbt0,sum0
      integer iphi,ibt,nphi,nbt,nbt1,nbt2,nbt0

      include 'pi.f'
      include 'ion.f'
      include 'beam.f'
      include 'p0Xn.f'
      include 'btmaxmin.f'
      include 'opacpbpars.f'
      include 'opacpbppars.f'


c      open(10,file='screen11a.dat')

c      nbt=10000
      if(beam.eq.'ion')nbt=ioppb
      if(beam.eq.'ionp')nbt=ioppbp

      nbt=3000
      nbt0=100

      if(ionbreakup)then
         if(fAA.eq.'00')then
            btmin=0d0
         elseif(fAA.eq.'0X'.or.fAA.eq.'X0')then
            btmax=1d6*rzg
            if(qt.gt.0.001d0)then
               btmax=1000d0*rzg/qt
            endif
            nbt=100000
         endif
      else
         btmin=0d0
      endif

      b0=4d0*rzg
      hbt0=(b0-btmin)/dble(nbt0)
      sum0=0d0

      do ibt=1,nbt0

         bt=btmin+(dble(ibt)-0.5d0)*hbt0

         if(beam.eq.'ion')wt=1d0-opacpbint(bt)
         if(beam.eq.'ionp')wt=1d0-opacpbpint(bt)

         wt=wt*besj0(qt*bt)
         wt=-wt*bt*hbt0/2d0/pi

         sum0=sum0+wt

      enddo

      btmin=b0

      hbt=(btmax-btmin)/dble(nbt)
      sum=0d0

      do ibt=1,nbt

         bt=btmin+(dble(ibt)-0.5d0)*hbt

         if(beam.eq.'ion')wt=1d0-opacpbint(bt)
         if(beam.eq.'ionp')wt=1d0-opacpbpint(bt)

c         if(beam.eq.'ion')wt=1d0-opacpbarr(2,ibt)
c         if(beam.eq.'ionp')wt=1d0-opacpbparr(2,ibt)

c$$$         if(bt.gt.2d0*rzg)then
c$$$            wt=0d0
c$$$         else
c$$$            wt=1d0
c$$$         endif

         wt=wt*besj0(qt*bt)
         wt=-wt*bt*hbt/2d0/pi

         sum=sum+wt

      enddo

      sum=sum+sum0

      screen=sum

      return
      end
