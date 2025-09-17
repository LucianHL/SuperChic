      function opacpb_inel(btt)
      implicit none
      double precision btmax,hphi,hbt,bt,phi,btp,btt
      integer iphi,ibt,nbt,nphi
      double precision sum,wt
      double precision btpx,btpy,btx,bty
      double precision opacpb_inel,opacpb,rhoxyint
      integer nbtx,nbty,ibty,ibtx
      double precision hbtx,hbty,opacpbint
      double precision breakup_prob

      include 'pi.f'
      include 'ion.f'
      include 'rho0.f'
      include 'ionqcd.f'
      include 'qcd.f'
      include 'p0Xn.f'
      

      btmax=1.5d0*rzg

      nphi=20
      nbt=20

      nbtx=20
      nbty=20

      hbtx=btmax/dble(nbtx)
      hbty=btmax/dble(nbty)

      hphi=2d0*pi/dble(nphi)
      hbt=btmax/dble(nbt)

      sum=0d0

c      do iphi=1,nphi
c      do ibt=1,nbt

c         bt=(dble(ibt)-0.5d0)*hbt
c         phi=(dble(iphi)-0.5d0)*hphi

c         btx=bt*dcos(phi)
c         bty=bt*dsin(phi)

      do ibtx=-nbtx,nbtx
      do ibty=-nbty,nbty

         btx=(dble(ibtx))*hbtx
         bty=(dble(ibty))*hbty
         bt=dsqrt(btx**2+bty**2)

         btpx=btx+btt
         btpy=bty

         btp=dsqrt(btpx**2+btpy**2)

c         if(opacpb(btp).gt.100d0)then
c         wt=0d0
c         else
c         wt=rhoxyint(1,bt)*dexp(-opacpb(btp))
c         endif
c         wt=rhoxyint(1,bt)*opacpbint(btp)**2

         if(qcd)then
         wt=(rhoxyint(1,bt)+rhoxyint(2,bt))*opacpbint(btp)**2
         else
         wt=rhoxyint(1,bt)*opacpbint(btp)**2
         endif
c         wt=wt*hbt*hphi*bt
         wt=wt*hbtx*hbty

         if(ionbreakup)wt=wt*breakup_prob(btp)

c         print*,bt/rzg,btp/rzg,wt,rhoxyint(1,bt)+rhoxyint(2,bt)
c     &   ,opacpbint(btp)

         sum=sum+wt

c         if(wt.gt.1d-10)print*,bt/rzg,btp,wt,dsqrt(sum/an)

      enddo
      enddo


      if(qcd)then
      opacpb_inel=dsqrt(sum/an)
      else
      opacpb_inel=dsqrt(sum/az)
      endif

      return
      end

      function opacpb_inel_incohqcd(btt)
      implicit none
      double precision btmax,hphi,hbt,bt,phi,btp,btt
      integer iphi,ibt,nbt,nphi
      double precision sum,wt
      double precision btpx,btpy,btx,bty,btx1,bty1,bt1
      double precision opacpb_inel_incohqcd,opacpb,rhoxyint
      integer nbtx,nbty,ibty,ibtx
      integer ibty1,ibtx1
      double precision hbtx,hbty,opacpbint
      double precision breakup_prob

      include 'pi.f'
      include 'ion.f'
      include 'rho0.f'
      include 'ionqcd.f'
      include 'qcd.f'
      include 'p0Xn.f'
      

      btmax=2d0*rzg

      nphi=10
      nbt=20

      nbtx=20
      nbty=20

      hbtx=btmax/dble(nbtx)
      hbty=btmax/dble(nbty)

      hphi=2d0*pi/dble(nphi)
      hbt=btmax/dble(nbt)

      sum=0d0

      do ibtx=-nbtx,nbtx
      do ibty=-nbty,nbty
      do ibtx1=-nbtx,nbtx
      do ibty1=-nbty,nbty

         btx=(dble(ibtx))*hbtx
         bty=(dble(ibty))*hbty

         btx1=(dble(ibtx1))*hbtx
         bty1=(dble(ibty1))*hbty

         bt=dsqrt(btx**2+bty**2)
         bt1=dsqrt(btx1**2+bty1**2)

         btpx=btx+btt+btx1
         btpy=bty+bty1

         btp=dsqrt(btpx**2+btpy**2)

         wt=(rhoxyint(1,bt)+rhoxyint(2,bt))*opacpbint(btp)**2
         wt=wt*(rhoxyint(1,bt1)+rhoxyint(2,bt1))


         wt=wt*hbtx*hbty*hbtx*hbty


         sum=sum+wt

c         if(dsqrt(wt/an**2).gt.1d-4)then
c         print*,rhoxyint(1,bt)+rhoxyint(2,bt)
c     &   ,rhoxyint(1,bt1)+rhoxyint(2,bt1),opacpbint(btp)
c         print*,bt/rzg,bt1/rzg,btp/rzg,dsqrt(wt/an**2),dsqrt(sum/an**2)
c         print*,''
c         endif

      enddo
      enddo
      enddo
      enddo


      if(qcd)then
      opacpb_inel_incohqcd=dsqrt(sum/an**2)
      else
      opacpb_inel_incohqcd=dsqrt(sum/az**2)
      endif

      return
      end


      function breakup_prob(bt)
      implicit none
      double precision bt,p0,pX,p1,pgdrint
      double precision breakup_prob

      include 'p0Xn.f'
      include 'ion.f'


c      if(fAA.eq.'00'.or.fAA.eq.'A0')then
c         if(bt.gt.rzg*1.5d0)then
c            p0=dexp(-pgdrint(3,dlog(bt)))
c         else
c            p0=0d0
c         endif
c      else
c         p0=dexp(-pgdrint(3,dlog(bt)))
c         pX=1d0-p0
c         p1=pgdrint(2,dlog(bt))*p0
c      endif

c      print*,'testa',bt
c      print*,bt,pgdrint(3,dlog(bt))
c      print*,'testb'


      if(bt.gt.rzg*1.5d0)then
         p0=dexp(-pgdrint(3,dlog(bt)))
      else
         p0=0d0
      endif
c      p0=dexp(-pgdrint(3,dlog(bt)))

      if(fAA.eq.'A0')breakup_prob=p0


      return
      end



      function opacpb(btt)
      implicit none
      double precision btmax
      double precision opacpb,wt,sum,rhoxy1,rhoxy2
      double precision phi1,phi2,hphi,hbt
      double precision btx,bty,bt2x,bt2y,bt1y,bt1x,bt2,bt1,bt,btt
      double precision rhoxyint,opacpint
      integer iphi1,iphi2,ibt1,ibt2
      integer nphi,nbt

      include 'pi.f'
      include 'ion.f'
      include 'rho0.f'
      include 'ionqcd.f'
      include 'qcd.f'
      include 'ion_inel.f'

      btmax=2d0*rzg
      nphi=20
      nbt=100

c      nphi=10
c      nbt=50

      hphi=2d0*pi/dble(nphi)
      hbt=btmax/dble(nbt)

      sum=0d0

      do 800 iphi1=1,nphi
      do 800 iphi2=1,nphi
      do 800 ibt1=1,nbt
      do 800 ibt2=1,nbt

         bt1=(dble(ibt1)-0.5d0)*hbt
         bt2=(dble(ibt2)-0.5d0)*hbt
         phi1=(dble(iphi1)-0.5d0)*hphi
         phi2=(dble(iphi2)-0.5d0)*hphi

         bt1x=bt1*dcos(phi1)
         bt1y=bt1*dsin(phi1)
         bt2x=bt2*dcos(phi2)
         bt2y=bt2*dsin(phi2)

         btx=btt-bt1x+bt2x
         bty=-bt1y+bt2y
         bt=dsqrt(btx**2+bty**2)

         if(opacpint(bt).gt.20d0)then
            wt=2d0
         else
            if(qcd)then
               if(ionqcd.eq.'coh')wt=2d0*(1d0-dexp(-opacpint(bt)/2d0))
               if(ionqcd.eq.'incoh')wt=(1d0-dexp(-opacpint(bt)))
            else
               wt=2d0*(1d0-dexp(-opacpint(bt)/2d0))
               if(ion_inel.and.ion_incoh_type.eq.'inel')then
               wt=(1d0-dexp(-opacpint(bt)))
               endif
            endif
         endif

c         rhoxy1=rhoxyint(1,bt1)+rhoxyint(2,bt1)
c         rhoxy2=rhoxyint(1,bt2)+rhoxyint(2,bt2)

         rhoxy1=rhoxyint(3,bt1)+rhoxyint(4,bt1)
         rhoxy2=rhoxyint(3,bt2)+rhoxyint(4,bt2)


         wt=wt*bt1*bt2*hphi**2*hbt**2
         wt=wt*rhoxy1*rhoxy2

         sum=sum+wt

 800  enddo

      opacpb=sum

      return
      end
