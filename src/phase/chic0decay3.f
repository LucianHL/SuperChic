ccc   spin correlations for chi_0 decay to mu+mu- + gamma via intermediate
ccc   vector meson
      subroutine chic0decay3(wtt)
      implicit none
      double precision mvec,q8q9,sdot,wt1,wtt
      complex*16 ze7q6,ze7q8,zce7q8,zce7q6,zce7e7
      complex*16 rho1psi(3,3)
      integer h,i,j,k
      complex*16 epsi1(3,4),epsi(4),cepsi(4)
      double precision q5(4),q6(4),q7(4),q8(4),q9(4)

      include 'polvecs.f'
      include 'mom.f'
      include 'quarkonia.f'

cccccccccccccccc

      do i=1,4
         q7(i)=q(i,7)
         q6(i)=q(i,6)
         q5(i)=q(i,5)
      enddo

      call genpol1(7,epsi1)
      mvec=dsqrt(q(4,7)**2-q(3,7)**2-q(2,7)**2-q(1,7)**2)
      
ccccccccccccccccccccc 

      

        do k=1,3
            do j=1,3
               
               do h=1,4
                  epsi(h)=epsi1(k,h)
                  cepsi(h)=conjg(epsi1(j,h))
               enddo
               
               call cdot(epsi,q6,ze7q6)
               call cdot(cepsi,q6,zce7q6)
               call ccdot(epsi,cepsi,zce7e7)

c     J/psi density matrix (unnormalised)
               
               rho1psi(k,j)=(-zce7e7*sdot(q6,q7)**2-mvec**2
     &              *ze7q6*zce7q6)/(2d0*sdot(q6,q7)**2)

            enddo
         enddo      

ccccccccccccccc

      do k=1,4
         q8(k)=q(k,8)
         q9(k)=q(k,9)
      enddo

      q8q9=sdot(q8,q9)

      wt1=0d0
    
      do j=1,3
         do k=1,3

            do h=1,4
               epsi(h)=epsi1(k,h)
               cepsi(h)=conjg(epsi1(j,h))
            enddo

               call cdot(epsi,q8,ze7q8)
               call cdot(cepsi,q8,zce7q8)
               call ccdot(epsi,cepsi,zce7e7)

            wt1=wt1+(-(q8q9+mmu**2)*zce7e7-2d0*ze7q8*
     &zce7q8)*rho1psi(k,j)/(2d0*(q8q9+2d0*mmu**2))*3d0

         enddo
      enddo

      wtt=wt1

      return
      end
