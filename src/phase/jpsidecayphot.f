ccc   spin correlations for j_psi decay to mu+mu- in photoproduction
ccc   (assuming SCHC)
      subroutine jpsidecayphot(wtt)
      implicit none
      double precision wtt,wt1,sh,q8q9
      complex*16 zq8epsi1,zq8cepsi1p,zepsi1cepsi1p
      complex*16 rho1psi(3,3)
      integer h,i,j,k
      complex*16 epsi1(3,4),epsi(4),cepsi(4)
      double precision q8(4)
      double precision pboo(4),plb(4),pcm(4)

      include 'polvecs.f'
      include 'mom.f'
      include 'quarkonia.f'
      include 'bpsi.f'
      include 'mres.f'
      include 'partonmom4.f'

cccc  psi density matrix:  3=long.
cccc                       1,2=transv.
cccc  default given by SCHC

      do k=1,3
         do j=1,3
            rho1psi(k,j)=0d0
         enddo
      enddo

      rho1psi(3,3)=0d0
      rho1psi(2,2)=(1d0-rho1psi(3,3))/2d0
      rho1psi(1,1)=rho1psi(2,2)

cccccccc
cccccccc    Boost to gamma-p rest frame (where SCHC holds)
cccccccc

      if(prot.eq.1)then
         do i=1,3
            pboo(i)=-q(i,1)+q(i,3)-q(i,2)
         enddo
         pboo(4)=q(4,1)-q(4,3)+q(4,2)
      elseif(prot.eq.2)then
         do i=1,3
            pboo(i)=-q(i,2)+q(i,4)-q(i,1)
         enddo
         pboo(4)=q(4,2)-q(4,4)+q(4,1)
      endif

      sh=dsqrt(pboo(4)**2-pboo(3)**2-pboo(2)**2-pboo(1)**2)

      do j=1,4
         pcm(j)=q(j,5)
      enddo
      call boost(sh,pboo,pcm,plb)  
      do j=1,4
         q(j,16)=plb(j)
      enddo

      do j=1,4
         pcm(j)=paa(j)
      enddo          
      do i=1,3
         pboo(i)=q(i,16)
      enddo
      pboo(4)=q(4,16)
      sh=dsqrt(pboo(4)**2-pboo(3)**2-pboo(2)**2-pboo(1)**2)

      call boost(sh,pboo,pcm,plb)  
      do j=1,4
         q(j,18)=plb(j)
      enddo

      call genpol1(16,epsi1)

ccccccccccccccccccc 

      do k=1,4
         q8(k)=q(k,18)
      enddo

      q8q9=(mres**2-2d0*mmu**2)/2d0

      wt1=0d0

      do j=1,3
         do k=1,3
            
            do h=1,4
               epsi(h)=epsi1(k,h)
               cepsi(h)=conjg(epsi1(j,h))           
            enddo

               call cdot(epsi,q8,zq8epsi1)
               call cdot(cepsi,q8,zq8cepsi1p)
               call ccdot(epsi,cepsi,zepsi1cepsi1p)

            wt1=wt1+dble((-(q8q9+mmu**2)*zepsi1cepsi1p-2d0*zq8epsi1*
     &zq8cepsi1p)*rho1psi(k,j)/(2d0*(q8q9+2d0*mmu**2))*3d0)

         enddo
      enddo

      wtt=wt1

      return
      end
