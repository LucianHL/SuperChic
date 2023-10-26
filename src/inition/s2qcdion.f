      subroutine s2qcdion
      implicit none
      double precision wt1,wt,sum,sum1,out
      double precision phi1,phi2,hbt,hphi
      double precision btx,bty,btmax,bt2x,bt2y,bt1y,bt1x,bt1,bt2,bt
      double precision tpqcdint,opacpbint
      integer iphi1,iphi2,ibt1,ibt2
      integer nphi,nbt

      include 'pi.f'
      include 'ion.f'
      include 'rho0.f'
      include 's2qcd.f'
      
      btmax=2d0*rzg
      nphi=20
      nbt=200

      hphi=2d0*pi/dble(nphi)
      hbt=btmax/dble(nbt)

      sum=0d0
      sum1=0d0
      
      do 800 iphi1=1,nphi
      do 800 ibt1=1,nbt
      do 800 ibt2=1,nbt

         bt1=(dble(ibt1)-0.5d0)*hbt
         bt2=(dble(ibt2)-0.5d0)*hbt
         phi1=(dble(iphi1)-0.5d0)*hphi
         phi2=(dble(iphi2)-0.5d0)*hphi

         phi2=0d0
         
         bt1x=bt1*dcos(phi1)
         bt1y=bt1*dsin(phi1)
         bt2x=bt2*dcos(phi2)
         bt2y=bt2*dsin(phi2)

         btx=bt1x-bt2x
         bty=bt1y-bt2y
         bt=dsqrt(btx**2+bty**2)

         wt=opacpbint(bt)**2
         wt=wt*bt1*bt2*hphi*hbt**2*2d0*pi
         wt=wt*(tpqcdint(bt1)*tpqcdint(bt2))

         wt1=bt1*bt2*hphi*hbt**2*2d0*pi
         wt1=wt1*(tpqcdint(bt1)*tpqcdint(bt2))
         
         sum=sum+wt
         sum1=sum1+wt1

 800  enddo
      
      out=sum/sum1
      
      s2qcd=out
      
      return
      end
