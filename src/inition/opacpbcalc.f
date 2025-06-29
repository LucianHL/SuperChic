      subroutine opacpbcalc
      implicit none
      double precision sigin,px,p1,p0,opac,sum1t,sum1ta
      double precision hb,bt,pgdrint,opacpb,prho,pll
     &   ,sum_cs,sum_cs1,sum_cs2,sum_cs3
      integer i
      logical readop

      include 'opacpbpars.f'
      include 'ion.f'
      include 'pi.f'
      include 'ionqcd.f'
      include 'qcd.f'
      include 'p0Xn.f'
      include 'btmaxmin.f'

      if(pAAvar)then
         if(ifaa.eq.1)then  ! inclusive
            ionbreakup=.false.
         endif
         if(ifaa.eq.2)then
            ionbreakup=.true.
            faa='00'
         endif
         if(ifaa.eq.3)then
            ionbreakup=.true.
            faa='XX'
         endif
      endif

      readop=.false.

      if(ionbreakup)then ! need bit more precision as wider b_t integration
         ioppb=400
         ioppb=2000
      else
         ioppb=200
      endif

      if(ionbreakup)then
         btmax=20d0*rzg
         btmax=200d0*rzg

cccccc

         btmin=rzg*1.6d0


         if(fAA.eq.'11')then
            btmax=rzg*200d0
            ioppb=3000
         elseif(fAA.eq.'00')then
            if(wrho)then
            btmax=rzg*200d0
            ioppb=3000
            else
            btmax=rzg*150d0
            btmin=rzg*1.7d0
            ioppb=1500
            endif
         elseif(fAA.eq.'A0')then
            if(wrho)then
            btmax=rzg*200d0
            ioppb=3000
            else
            btmax=rzg*150d0
            btmin=rzg*1.7d0
            ioppb=1500
            endif
         elseif(fAA.eq.'AA')then
            if(wrho)then
            btmax=rzg*200d0
            ioppb=3000
            else
            btmax=rzg*4d0
            btmin=rzg*1.5d0
            ioppb=25
            endif
         elseif(fAA.eq.'XX')then
            btmax=rzg*200d0
            ioppb=3000
         elseif(fAA.eq.'1X'.or.fAA.eq.'X1')then
            btmax=rzg*200d0
            ioppb=3000
         elseif(fAA.eq.'01')then
            btmax=rzg*200d0
            ioppb=3000
         elseif(fAA.eq.'A1')then
            btmax=rzg*200d0
            ioppb=3000
         elseif(fAA.eq.'AX')then
            btmax=rzg*200d0
            ioppb=3000
         endif



      else
         btmax=rzg*4d0
         btmin=rzg*1.5d0
         ioppb=25
      endif

      hb=(btmax-btmin)/dble(ioppb)

      if(ionqcd.eq.'incoh')sigin=6.5d0 ! fm^2
      if(ionqcd.eq.'coh')sigin=6.74d0 ! fm^2

      do i=1,ioppb

         bt=btmin+(dble(i)-0.5d0)*hb

         opacpbarr(1,i)=bt

         if(bt.lt.rzg*1.5d0)then
            opac=200d0
         elseif(bt.lt.4d0*rzg)then ! above this no point evaluating, ~ 0
            opac=opacpb(bt)
         else
            opac=0d0
         endif


         if(qcd)opac=opac*(1d0-1d0/nshell)**2
         if(opac.gt.100d0)then
            opacpbarr(2,i)=0d0
         else
            opacpbarr(2,i)=dexp(-opac/2d0)
         endif

         if(ionbreakup)then


            if(wrho)then
               p0=dexp(-pgdrint(3,dlog(bt)))
               pX=1d0-p0
               p1=pgdrint(2,dlog(bt))*p0
            elseif(fAA.eq.'00'.or.fAA.eq.'A0')then
               if(bt.gt.rzg*1.5d0)then
                  p0=dexp(-pgdrint(3,dlog(bt)))
               else
                  p0=0d0
               endif
            else
               p0=dexp(-pgdrint(3,dlog(bt)))
               pX=1d0-p0
               p1=pgdrint(2,dlog(bt))*p0
            endif

            if(wrho)then
               prho=pgdrint(4,dlog(bt))
               pll=pgdrint(5,dlog(bt))
            endif

            approx=.false.


            if(wrho)then

            sum1t=pgdrint(4,dlog(btmax))/bt**2*btmax**2

            else

            sum1t=0.75d0*az*(an-az)/an**(2d0/3d0)
     &              /0.389389d0*az**2/pi**2/bt**2/137d0

            endif

            if(approx)then
               sum1t=0.75d0*az*(an-az)/an**(2d0/3d0)
     &              /0.389389d0*az**2/pi**2/bt**2/137d0
               p0=dexp(-sum1t)
               pX=1d0-p0
               p1=sum1t*p0
            endif

            if(fAA.eq.'00')then
               if(wrho)then
               opacpbarr(3,i)=1d0-opacpbarr(2,i)*dsqrt(sum1t)
               opacpbarr(2,i)=1d0-opacpbarr(2,i)*(dsqrt(prho)*p0
     &              -dsqrt(sum1t))
               else
               opacpbarr(2,i)=opacpbarr(2,i)*p0
               endif
            elseif(fAA.eq.'XX')then
               if(wrho)then
               opacpbarr(2,i)=1d0-opacpbarr(2,i)*pX*dsqrt(prho)
               else
               opacpbarr(2,i)=1d0-opacpbarr(2,i)*pX
               endif
            elseif(fAA.eq.'11')then
               if(wrho)then
               opacpbarr(2,i)=1d0-opacpbarr(2,i)*p1*dsqrt(prho)
               else
               opacpbarr(2,i)=1d0-opacpbarr(2,i)*p1
               endif
            elseif(fAA.eq.'0X')then
               if(wrho)then
               opacpbarr(2,i)=1d0-opacpbarr(2,i)*dsqrt(pX*p0*prho)
               else
               opacpbarr(2,i)=1d0-opacpbarr(2,i)*dsqrt(pX*p0)
               endif
            elseif(fAA.eq.'01')then
               if(wrho)then
               opacpbarr(2,i)=1d0-opacpbarr(2,i)*dsqrt(p1*p0*prho)
               else
               opacpbarr(3,i)=1d0-opacpbarr(2,i)*dsqrt(sum1t)
               opacpbarr(2,i)=1d0-opacpbarr(2,i)*(dsqrt(p1*p0)
     &              -dsqrt(sum1t))
               endif
            elseif(fAA.eq.'1X')then
               if(wrho)then
               opacpbarr(2,i)=1d0-opacpbarr(2,i)*dsqrt(p1*pX*prho)
               else
               opacpbarr(2,i)=1d0-opacpbarr(2,i)*dsqrt(p1*pX)
               endif
            elseif(fAA.eq.'AA')then
               if(wrho)then
               opacpbarr(3,i)=1d0-opacpbarr(2,i)*dsqrt(sum1t)
               opacpbarr(2,i)=1d0-opacpbarr(2,i)*(dsqrt(prho)
     &              -dsqrt(sum1t))
c               print*,bt,opacpbarr(3,i),opacpbarr(2,i)
               else
               opacpbarr(2,i)=opacpbarr(2,i)
               endif
            elseif(fAA.eq.'A0')then
               if(wrho)then
               opacpbarr(3,i)=1d0-opacpbarr(2,i)*dsqrt(sum1t)
               opacpbarr(2,i)=1d0-opacpbarr(2,i)*(dsqrt(prho*p0)
     &              -dsqrt(sum1t))
               else
               opacpbarr(2,i)=opacpbarr(2,i)*dsqrt(p0)
               endif
            elseif(fAA.eq.'AX')then
               if(wrho)then
               opacpbarr(2,i)=1d0-opacpbarr(2,i)*dsqrt(pX*prho)
               else
               opacpbarr(3,i)=1d0-opacpbarr(2,i)*dsqrt(sum1t)
               opacpbarr(2,i)=1d0-opacpbarr(2,i)*(dsqrt(pX)
     &              -dsqrt(sum1t))
               endif
            elseif(fAA.eq.'A1')then
               if(wrho)then
               opacpbarr(2,i)=1d0-opacpbarr(2,i)*dsqrt(p1*prho)
               else
               opacpbarr(3,i)=1d0-opacpbarr(2,i)*dsqrt(sum1t)
               opacpbarr(2,i)=1d0-opacpbarr(2,i)*(dsqrt(p1)
     &              -dsqrt(sum1t))
               endif
            else
               print*,'fAA option not allowed - STOP'
               STOP 1
            endif

         endif

      enddo


      return
      end
