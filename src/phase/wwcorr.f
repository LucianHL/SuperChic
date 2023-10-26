ccc   spin correlations in W+W- production (leptonic decays)
      subroutine wwcorr(wt,wtt)
      implicit none
      double precision sint2,sint1,rhowpp,rhowpm,rhowp0,rhowmp,rhowm0
     &     ,rhowmm,cost1,cost2,wtt
      double precision rhow(9)
      integer mm
      complex*16 wt(10),rhoww(9)
      logical rf
      double precision rhowp_rf(3),rhowm_rf(3)

      include 'ewpars.f'
      include 'pi.f'
      include 'partonmom2.f'
      include 'partonmom4.f'
      include 'mom.f'

      do mm=1,9
         rhoww(mm)=cdabs(wt(mm))**2
      enddo

       rf=.true.

      cost1=(q(1,6)*paa(1)+q(2,6)*paa(2)+q(3,6)*paa(3))
     &/dsqrt((paa(1)**2+paa(2)**2+paa(3)**2)*
     &     (q(1,6)**2+q(2,6)**2+q(3,6)**2))
      sint1=dsqrt(1d0-cost1**2)
      
      cost2=(q(1,7)*pbb(1)+q(2,7)*pbb(2)+q(3,7)*pbb(3))
     &/dsqrt((pbb(1)**2+pbb(2)**2+pbb(3)**2)*
     &     (q(1,7)**2+q(2,7)**2+q(3,7)**2))
      sint2=dsqrt(1d0-cost2**2)
      
      rhowp0=-sint1*dsqrt(3d0/2d0)
      rhowpm=(1d0+cost1)*dsqrt(3d0/4d0)
      rhowpp=(1d0-cost1)*dsqrt(3d0/4d0)
      rhowm0=-sint2*dsqrt(3d0/2d0)
      rhowmp=(1d0+cost2)*dsqrt(3d0/4d0)
      rhowmm=(1d0-cost2)*dsqrt(3d0/4d0)
      
      rhow(1)=rhowpp*rhowmp
      rhow(2)=rhowpp*rhowmm
      rhow(3)=rhowpm*rhowmp
      rhow(4)=rhowpm*rhowmm
      rhow(5)=rhowp0*rhowmp  
      rhow(6)=rhowp0*rhowmm      
      rhow(7)=rhowpp*rhowm0
      rhow(8)=rhowpm*rhowm0
      rhow(9)=rhowp0*rhowm0  

      if(rf)then
      
         call dmat_rf(1,rhowp_rf)
         call dmat_rf(2,rhowm_rf)

         rhow(1)=dsqrt(rhowp_rf(1)*rhowm_rf(1))
         rhow(2)=dsqrt(rhowp_rf(1)*rhowm_rf(2))
         rhow(3)=dsqrt(rhowp_rf(2)*rhowm_rf(1))
         rhow(4)=dsqrt(rhowp_rf(2)*rhowm_rf(2))
         rhow(5)=dsqrt(rhowp_rf(3)*rhowm_rf(1))
         rhow(6)=dsqrt(rhowp_rf(3)*rhowm_rf(2))
         rhow(7)=dsqrt(rhowp_rf(1)*rhowm_rf(3))
         rhow(8)=dsqrt(rhowp_rf(2)*rhowm_rf(3))       
         rhow(9)=dsqrt(rhowp_rf(3)*rhowm_rf(3))

      endif
      
      wtt=0d0

      do mm=1,9
         wtt=wtt+rhoww(mm)*rhow(mm)**2
      enddo

      return
      end

      subroutine dmat_rf_old(iw,rho_rf)
      implicit none
      double precision mw
      double precision rho_rf(3)
      integer i,j,mu,nu,beta,iw
      complex*16 ep(4),pl_e,pn_e,eps,zt,ep_epc
      
      include 'polwrf.f'
      include 'gmatrices.f'
      include 'partonmom4.f'
      include 'partonmom2.f'
      include 'zi.f'

      mw=dsqrt(p1(4)**2-p1(3)**2-p1(2)**2-p1(1)**2)

      
      do i=1,3

         if(i.eq.1)then
            do j=1,4
               ep(j)=echirf(iw,1,j)
            enddo
         elseif(i.eq.2)then
            do j=1,4
               ep(j)=echirf(iw,2,j)
            enddo
         else
            do j=1,4
               ep(j)=echirf(iw,3,j)
            enddo
         endif

         ep_epc=ep(4)*dconjg(ep(4))-ep(3)*dconjg(ep(3))
     &        -ep(2)*dconjg(ep(2))-ep(1)*dconjg(ep(1))

 
         
         if(iw.eq.1)then
            pl_e=paa(4)*ep(4)-paa(3)*ep(3)-paa(2)*ep(2)-paa(1)*ep(1)
            pn_e=mw*ep(4)-pl_e
         else
            pn_e=pbb(4)*ep(4)-pbb(3)*ep(3)-pbb(2)*ep(2)-pbb(1)*ep(1)
            pl_e=mw*ep(4)-pn_e
         endif

         eps=0d0
         
         do mu=1,4
         do nu=1,4
         do beta=1,4

            if(iw.eq.1)then
               zt=zi*e_(mu,nu,4,beta)*paa(mu)*ep(nu)*mw*
     &              dconjg(ep(beta))
            else
               zt=zi*e_(4,nu,mu,beta)*pbb(mu)*ep(nu)*mw*
     &              dconjg(ep(beta))
            endif

            if(mu.lt.4)zt=-zt
            if(nu.lt.4)zt=-zt
            if(beta.lt.4)zt=-zt

            eps=eps+zt
            
         enddo
         enddo
         enddo

         rho_rf(i)=4d0*(2d0*dreal(pl_e*dconjg(pn_e))-ep_epc*mw**2/2d0)
     &        -4d0*eps
         rho_rf(i)=rho_rf(i)/mw**2*3d0/4d0

      enddo
      
      return
      end

      
      subroutine dmat_rf(iw,rho_rf)
      implicit none
      double precision ml,mw
      double precision rho_rf(3)
      integer i,j,mu,nu,beta,iw
      complex*16 ep(4),pf_e,paf_e,eps,zt,ep_epc,pw_e
      
      include 'polwrf.f'
      include 'gmatrices.f'
      include 'partonmom4.f'
      include 'partonmom2.f'
      include 'zi.f'
      include 'mom.f'
      include 'wwpars.f'

      mw=dsqrt(p1(4)**2-p1(3)**2-p1(2)**2-p1(1)**2)
      
      do i=1,3

         if(i.eq.1)then
            do j=1,4
               ep(j)=echirf(iw,1,j)
            enddo
         elseif(i.eq.2)then
            do j=1,4
               ep(j)=echirf(iw,2,j)
            enddo
         else
            do j=1,4
               ep(j)=echirf(iw,3,j)
            enddo
         endif

         ep_epc=ep(4)*dconjg(ep(4))-ep(3)*dconjg(ep(3))
     &        -ep(2)*dconjg(ep(2))-ep(1)*dconjg(ep(1))

         if(iw.eq.1)then ! w+
            pf_e=paa(4)*ep(4)-paa(3)*ep(3)-paa(2)*ep(2)-paa(1)*ep(1) ! fermion = neutrino
            paf_e=mw*ep(4)-pf_e ! antifermion = lepton (e,mu)
            ml=dsqrt(dabs(q(4,9)**2-q(3,9)**2-q(2,9)**2-q(1,9)**2)) ! more stable
         else ! w-
            paf_e=pbb(4)*ep(4)-pbb(3)*ep(3)-pbb(2)*ep(2)-pbb(1)*ep(1) ! antifermion = neutrino
            pf_e=mw*ep(4)-paf_e ! fermion = lepton (e,mu)
            ml=dsqrt(dabs(q(4,11)**2-q(3,11)**2-q(2,11)**2-q(1,11)**2)) ! more stable
         endif

         pw_e=mw*ep(4)
         
         eps=0d0
         
         do mu=1,4
         do nu=1,4
         do beta=1,4

            if(iw.eq.1)then
               zt=zi*e_(mu,nu,4,beta)*paa(mu)*ep(nu)*mw*
     &              dconjg(ep(beta))
            else
               zt=zi*e_(4,nu,mu,beta)*pbb(mu)*ep(nu)*mw*
     &              dconjg(ep(beta))
            endif

            if(mu.lt.4)zt=-zt
            if(nu.lt.4)zt=-zt
            if(beta.lt.4)zt=-zt

            eps=eps+zt
            
         enddo
         enddo
         enddo

         if(wgauge.eq.'unitary')then
         
         rho_rf(i)=4d0*(2d0*dreal(pf_e*dconjg(paf_e))
     &        +ep_epc*(ml**2-mw**2)/2d0)-4d0*eps
         rho_rf(i)=rho_rf(i)*3d0/4d0/(mw**2-ml**2/2d0-ml**4/mw**2/2d0)
    
         rho_rf(i)=dabs(rho_rf(i))

         else
         
         if(iw.eq.1)then
         
            rho_rf(i)=4d0*(2d0*dreal(pf_e*dconjg(paf_e))
     &           -2d0*dreal(pf_e*dconjg(pw_e))*ml**2/mw**2
     &           +pw_e*dconjg(pw_e)*ml**2/mw**2*(1d0-ml**2/mw**2)/2d0
     &           +ep_epc*(ml**2-mw**2)/2d0)-4d0*eps
           rho_rf(i)=rho_rf(i)*3d0/4d0/(mw**2-ml**2/2d0-ml**4/mw**2/2d0)

         else

            rho_rf(i)=4d0*(2d0*dreal(pf_e*dconjg(paf_e))
     &           +2d0*dreal(paf_e*dconjg(pw_e))*ml**2/mw**2
     &           +pw_e*dconjg(pw_e)*ml**2/mw**2*(1d0-ml**2/mw**2)/2d0
     &           +ep_epc*(ml**2-mw**2)/2d0)-4d0*eps
           rho_rf(i)=rho_rf(i)*3d0/4d0/(mw**2-ml**2/2d0-ml**4/mw**2/2d0)
            
         endif

         rho_rf(i)=dabs(rho_rf(i))

         endif
         
      enddo
      
      return
      end
