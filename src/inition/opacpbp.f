      function opacpbp(btt)
      implicit double precision(a-y)
      integer iphi1,iphi2,ibt1,ibt2
      integer nphi,nbt

      include 'pi.f'
      include 'ion.f'
      include 'rho0.f'
      include 'ionqcd.f'
      include 'qcd.f'
      
      rhop=0.1973d0**2/pi ! proton transverse area
      
      btmax=5d0*rzg
      nphi=20
      nbt=200
      
      hphi=2d0*pi/dble(nphi)
      hbt=btmax/dble(nbt)

      sum=0d0
      
      do 800 iphi1=1,nphi
      do 800 ibt1=1,nbt

         bt1=(dble(ibt1)-0.5d0)*hbt
         phi1=(dble(iphi1)-0.5d0)*hphi

         bt1x=bt1*dcos(phi1)
         bt1y=bt1*dsin(phi1)

         btx=btt-bt1x
         bty=-bt1y
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

         wt=wt*bt1*hphi*hbt
c         wt=wt*(rhoxyint(1,bt1)+rhoxyint(2,bt1))
         wt=wt*(rhoxyint(3,bt1)+rhoxyint(4,bt1))
            
         sum=sum+wt
         
 800  enddo
      
      opacpbp=sum
      
      return
      end
