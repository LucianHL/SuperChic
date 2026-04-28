      subroutine opacpbcalc
      implicit none
      double precision sigin,px,p1,p0,opac,sum1t,sum1ta,p0h
      double precision hb,bt,pgdrint,opacpb,prho,pll,opacpb1
     &   ,sum_cs,sum_cs1,sum_cs2,sum_cs3
      integer i,ioppba
      logical readop
      double precision opacpb_inel,opacpbint,opacpb_inel_incohqcd
      double precision test,test1

      include 'opacpbpars.f'
      include 'ion.f'
      include 'pi.f'
      include 'ionqcd.f'
      include 'qcd.f'
      include 'p0Xn.f'
      include 'btmaxmin.f'
      include 'ion_inel.f'
      include 'veto.f'


      if(pAAvar)then
         if(ifaa.eq.1)then  ! inclusive
c            ionbreakup=.false.
            ionbreakup=.true. ! NEW
            faa='AA'
            faa='AX'
            int_01=.true.
         endif
         if(ifaa.eq.2)then
            ionbreakup=.true.
            faa='00'
            faa='XX'
            int_01=.false.
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
c         ioppb=200
         ioppb=2000
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
            btmin=rzg*1.5d0
            ioppb=3000
            btmax=rzg*300d0
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
         elseif(fAA.eq.'0X')then
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
c         btmax=rzg*4d0
c         btmin=rzg*1.5d0
c         ioppb=25

         btmin=rzg*1.5d0
         ioppb=3000
         btmax=rzg*300d0

      endif

      if(ion_inel)then
         btmax=rzg*4d0
         btmin=rzg*1.5d0
         ioppb=25
      

      endif

      hb=(btmax-btmin)/dble(ioppb)

      if(ionqcd.eq.'incoh')sigin=6.5d0 ! fm^2
      if(ionqcd.eq.'coh')sigin=6.74d0 ! fm^2

      if(ion_inel)goto 444
      if(qcd.and.ionqcd.eq.'incoh')goto 444

      sum_cs=0d0

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

               p0=dexp(-pgdrint(3,dlog(bt)))
               pX=1d0-p0

            else
               p0=dexp(-pgdrint(3,dlog(bt)))
               pX=1d0-p0

               p1=pgdrint(2,dlog(bt))*p0
            endif

            if(veto)then
               p0h=dexp(-pgdrint(4,dlog(bt)))
            else
               p0h=1d0
            endif


            if(wrho)then
               prho=pgdrint(5,dlog(bt))
c               pll=pgdrint(5,dlog(bt))
            endif

            approx=.false.


            if(wrho)then

c            sum1t=pgdrint(4,dlog(btmax))/bt**2*btmax**2
            sum1t=pgdrint(5,dlog(btmax))/bt**2*btmax**2

            else

            sum1t=0.75d0*az*(an-az)/an**(2d0/3d0)
     &              /0.389389d0*az**2/pi**2/bt**2/137d0

            sum1t=pgdrint(3,dlog(btmax))*btmax**2/bt**2


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
               opacpbarr(2,i)=1d0-opacpbarr(2,i)*(dsqrt(prho)*p0*p0h  
     &              -dsqrt(sum1t))
               else
               opacpbarr(2,i)=opacpbarr(2,i)*p0*p0h  ! NEW
               endif
            elseif(fAA.eq.'XX')then
               if(wrho)then
               opacpbarr(2,i)=1d0-opacpbarr(2,i)*pX*dsqrt(prho)*p0h  
               else
               opacpbarr(2,i)=1d0-opacpbarr(2,i)*pX*p0h  ! NEW
               endif
            elseif(fAA.eq.'11')then
               if(wrho)then
               opacpbarr(2,i)=1d0-opacpbarr(2,i)*p1*dsqrt(prho)*p0h  
               else
               opacpbarr(2,i)=1d0-opacpbarr(2,i)*p1*p0h  
               endif
            elseif(fAA.eq.'0X')then
               if(wrho)then
               opacpbarr(2,i)=1d0-opacpbarr(2,i)*dsqrt(pX*p0*prho)*p0h  
               else
                opacpbarr(3,i)=1d0-opacpbarr(2,i)*dsqrt(sum1t)
                opacpbarr(2,i)=1d0-opacpbarr(2,i)*(dsqrt(pX*p0)*p0h
     &              -dsqrt(sum1t))
               endif
            elseif(fAA.eq.'01')then
               if(wrho)then
               opacpbarr(2,i)=1d0-opacpbarr(2,i)*dsqrt(p1*p0*prho)*p0h  
               else
               opacpbarr(3,i)=1d0-opacpbarr(2,i)*dsqrt(sum1t)
               opacpbarr(2,i)=1d0-opacpbarr(2,i)*(dsqrt(p1*p0)*p0h  
     &              -dsqrt(sum1t))
               endif
            elseif(fAA.eq.'1X')then
               if(wrho)then
               opacpbarr(2,i)=1d0-opacpbarr(2,i)*dsqrt(p1*pX*prho)*p0h  
               else
               opacpbarr(2,i)=1d0-opacpbarr(2,i)*dsqrt(p1*pX)*p0h  
               endif
            elseif(fAA.eq.'AA')then
               if(wrho)then
               opacpbarr(3,i)=1d0-opacpbarr(2,i)*dsqrt(sum1t)
               opacpbarr(2,i)=1d0-opacpbarr(2,i)*(dsqrt(prho)*p0h  
     &              -dsqrt(sum1t))
               else
                opacpbarr(2,i)=opacpbarr(2,i)*p0h  ! NEW
               endif
            elseif(fAA.eq.'A0')then
               if(wrho)then
               opacpbarr(3,i)=1d0-opacpbarr(2,i)*dsqrt(sum1t)
               opacpbarr(2,i)=1d0-opacpbarr(2,i)*(dsqrt(prho*p0)*p0h  
     &              -dsqrt(sum1t))
               else
               opacpbarr(2,i)=opacpbarr(2,i)*dsqrt(p0)*p0h
               endif
            elseif(fAA.eq.'AX')then
               if(wrho)then
               opacpbarr(2,i)=1d0-opacpbarr(2,i)*dsqrt(pX*prho)*p0h  
               else
               opacpbarr(3,i)=1d0-opacpbarr(2,i)*dsqrt(sum1t)
               opacpbarr(2,i)=1d0-opacpbarr(2,i)*(dsqrt(pX)*p0h
     &              -dsqrt(sum1t))
               endif
            elseif(fAA.eq.'A1')then
               if(wrho)then
               opacpbarr(2,i)=1d0-opacpbarr(2,i)*dsqrt(p1*prho)*p0h  
               else
               opacpbarr(3,i)=1d0-opacpbarr(2,i)*dsqrt(sum1t)
               opacpbarr(2,i)=1d0-opacpbarr(2,i)*(dsqrt(p1)*p0h  
     &              -dsqrt(sum1t))
               endif
            else
               print*,'fAA option not allowed - STOP'
               STOP 1
            endif

         endif

      enddo

      return



  444 continue


      if(ionbreakup)then

      if(fAA.eq.'A0')then
         btmax=rzg*150d0
         btmin=rzg*1.5d0
         ioppb=1500
      endif

      endif

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

         if(opac.gt.100d0)then
            opacpbarr(2,i)=0d0
         else
            opacpbarr(2,i)=dexp(-opac/2d0)
         endif


      enddo

      btmin=0d0
      hb=(btmax-btmin)/dble(ioppb)



      do i=1,ioppb

         bt=btmin+(dble(i)-0.5d0)*hb

         if(qcd)then
         if(ionqcd.eq.'incoh')then
         opacpbarr_temp(i)=opacpb_inel_incohqcd(bt)
         else
         opacpbarr_temp(i)=opacpb_inel(bt)
         endif
         else
         opacpbarr_temp(i)=opacpb_inel(bt)
         endif


         if(opacpbarr_temp(i).gt.1d0)opacpbarr_temp(i)=1d0

         if(ionbreakup)then

            if(fAA.eq.'00'.or.fAA.eq.'A0')then
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


         endif



      enddo


cccc  And now overewrite (standard opacity needed above)


      do i=1,ioppb
         bt=btmin+(dble(i)-0.5d0)*hb
         opacpbarr(1,i)=bt
         opacpbarr(2,i)=opacpbarr_temp(i)

      enddo


      return
      end
