      subroutine screencalc
      implicit none
      double precision qmin,qmax,qt,lqmin,lqmax,lq
      double precision screen,screen_01,screeninv
      double precision btmin,btmax,bt,hbt,opacpbint,output
      double precision p0,pX,pgdrint,orig,opac,opacpb
      integer i,itotb

      include 'scionpars.f'
      include 'p0Xn.f'
      include 'ion.f'
      include 'beam.f'

      qmin=1d-5
      qmax=1.5d0

      qmin=1d-7

      lqmin=dlog(qmin)
      lqmax=dlog(qmax)

      itot=400
      itot=800  ! NEW

      do i=1,itot+1

         lq=lqmin+(lqmax-lqmin)*dble(i-1)/dble(itot)
         qt=dexp(lq)

         scionarr(1,i)=lq

         if(pAAvar)then
c            scionarr(ifaa+1,i)=screen(qt)
            if(int_01)then
            scionarr(ifaa+1,i)=screen_01(qt)
            scionarr(ifaa+1,i)=scionarr(ifaa+1,i)+screen(qt)
            else
            scionarr(ifaa+1,i)=screen(qt)
            endif
         elseif(int_01)then
            scionarr(2,i)=screen_01(qt)
            scionarr(2,i)=scionarr(2,i)+screen(qt)
         else
            scionarr(2,i)=screen(qt)
         endif

      enddo

cccccccccc

      itotb=200

      btmin=0d0
      btmax=50d0*rzg
      hbt=(btmax-btmin)/dble(itotb)

      goto 777

      do i=1,itotb+1

         bt=btmin+hbt*(dble(i)-0.5d0)

         p0=dexp(-pgdrint(3,dlog(bt)))
         pX=1d0-p0
         opac=opacpb(bt)

         if(opac.gt.100d0)then
            opac=0d0
         else
            opac=dexp(-opac/2d0)
         endif

         orig=1d0-opac
         orig=opac
c         *dsqrt(p0)
         orig=1d0-opac*dsqrt(pX)
      

      enddo

777   continue

      return
      end

      subroutine screen_int(output)
      implicit none
      integer i,itot
      double precision screeningionint,output
      double precision bt,sum,kt,ktmin,ktmax,hkt,wt

      itot=800

      sum=0d0

      ktmin=0d0
      ktmax=1.5d0

      hkt=(ktmax-ktmin)/dble(itot)


      do i=1,itot

      kt=ktmin+hkt*(dble(i)-0.5d0)

c      print*,kt

      wt=kt*hkt
      wt=wt*kt
      wt=wt*screeningionint(kt)

      sum=sum+wt

      enddo

      output=sum**2

      return
      end




      function screeninv(bt)
      implicit none
      integer i,itot
      double precision screeninv,screeningionint
      double precision bt,sum,kt,ktmin,ktmax,hkt,wt

      include 'pi.f'

      itot=16000

      sum=0d0

      ktmin=0d0
      ktmax=1.5d0

      hkt=(ktmax-ktmin)/dble(itot)


      do i=1,itot

      kt=ktmin+hkt*(dble(i)-0.5d0)

c      print*,kt

      wt=2d0*pi*kt*hkt
      wt=wt*BESSEL_J0(kt*bt)
      wt=wt*screeningionint(kt)

      sum=sum+wt

      enddo

      screeninv=1d0+sum

      return
      end