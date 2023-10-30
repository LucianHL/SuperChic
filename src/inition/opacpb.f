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

      btmax=2d0*rzg
      nphi=20
      nbt=100

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
