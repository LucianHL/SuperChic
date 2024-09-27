      function tpqcd(rxy)
      implicit none
      double precision rzmax,hz,tbtmin,tbtmax,tpqcd
      double precision b0,rxy,rhoxyint,htbt
      integer itot
      integer nbt,nphi,nk

      include 'ion.f'
      include 'pi.f'
      include 'rho0.f'

      itot=10000

      rzmax=5d0*rzg
      hz=rzmax/dble(itot)

      nk=200

      nbt=200
      nphi=20

      b0=2d0

      rzmax=10d0*dsqrt(b0)

      tbtmin=b0*dexp(-rzmax**2/b0)
      tbtmax=b0
      htbt=(tbtmax-tbtmin)/dble(nbt)

c$$$      hbt=rzmax/dble(nbt)
c$$$      hphi=2d0*pi/dble(nphi)
c$$$
c$$$      sum=0d0
c$$$      sumn=0d0
c$$$
c$$$      ktmax=2d0
c$$$      hk=ktmax/dble(nk)
c$$$
c$$$      do ibt=1,nbt
c$$$         do iphi=1,nphi
c$$$
c$$$            phi=hphi*(dble(iphi)-0.5d0)
c$$$            bt=hbt*(dble(ibt)-0.5d0)
c$$$
c$$$            btx=bt*dcos(phi)
c$$$            bty=bt*dsin(phi)
c$$$
c$$$            btdif=(rxy-btx)**2+bty**2
c$$$
c$$$            wt=rhoxyint(1,dsqrt(btdif))+rhoxyint(2,dsqrt(btdif))
c$$$            wt=wt*dexp(-bt**2/b0/4d0)
c$$$            wt=wt/(4d0*pi*b0)
c$$$            wt=wt*bt*hbt*hphi
c$$$
c$$$            sum=sum+wt
c$$$
c$$$         enddo
c$$$      enddo
c$$$
c$$$      tpqcd=sum

      tpqcd=rhoxyint(1,rxy)+rhoxyint(2,rxy)

      return
      end
