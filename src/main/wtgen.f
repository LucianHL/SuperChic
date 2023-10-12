ccc   calls subprocess amplitude
      subroutine wtgen
      implicit none
      integer p
      complex*16 ztest,zpp,zmm,zpm,zmp
      double precision mu

      include 'polarization.f'
      include 'zarr.f'
      include 'mandelstam.f'
      include 'vars.f'
      include 'proc.f'
      include 'mq.f'

      do p=1,pol

      if(proc.eq.1)then
            call higgs(mx,zpp,zmm,zpm,zmp)
         elseif(proc.eq.2)then
            call setmu(mu)    
            call gamgam(p,mu,uh,th,zpp,zmm,zpm,zmp)
         elseif(proc.eq.3)then
            call gg(p,mx,uh,th,zpp,zmm,zpm,zmp)
         elseif(proc.eq.4.or.proc.eq.5.or.proc.eq.6)then
            call qq(p,mx,mq,uh,th,zpp,zmm,zpm,zmp)
         elseif(proc.eq.7)then
            call ggg(p,mx,zpp,zmm,zpm,zmp)
         elseif(proc.eq.8)then
            call qqg(p,mx,zpp,zmm,zpm,zmp)
         elseif(proc.eq.9.or.proc.eq.10.or.proc.eq.11.or.proc.eq.12)then
            call pipi(p,mx,th,zpp,zmm,zpm,zmp)
         elseif(proc.eq.13)then 
            call rhorho(p,mx,th,zpp,zmm,zpm,zmp) 
         elseif(proc.eq.14)then 
            call etaeta(p,mx,uh,th,zpp,zmm,zpm,zmp) 
         elseif(proc.eq.15)then 
            call etaetap(p,mx,uh,th,zpp,zmm,zpm,zmp) 
         elseif(proc.eq.16)then 
            call etapetap(p,mx,uh,th,zpp,zmm,zpm,zmp) 
         elseif(proc.eq.17)then
            call rhorho(p,mx,th,zpp,zmm,zpm,zmp) 
          elseif(proc.eq.18)then
             call setmu(mu)     
             call psipsi(p,mu,mx,uh,th,zpp,zmm,zpm,zmp)
         elseif(proc.eq.19)then
             call setmu(mu)    
             call psipsip(p,mu,mx,uh,th,zpp,zmm,zpm,zmp)
         elseif(proc.eq.20)then
             call setmu(mu)     
             call psippsip(p,mu,mx,uh,th,zpp,zmm,zpm,zmp)
         endif

         zarr(1,p)=zpp
         zarr(2,p)=zmm
         zarr(3,p)=zpm
         zarr(4,p)=zmp

         enddo

      return
      end
