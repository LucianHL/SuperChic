ccc   spin correlations for chi_1 decay to mu+mu- + gamma via intermediate
ccc   vector meson
      subroutine chic1decay3(wt,wtt)
      implicit none
      double precision wtt,sdot,q8q9,wt1,pnorm,mvec
      complex*16 ze7q8,ze7q6,ze5q6,ze5e7,zce7q8,zce7q6,zce7e7,zce5q6,
     &     zce5e5,zce5ce7
      complex*16 wt(10)
      complex*16 rho1psi(3,3),rho1chi(3,3)
      integer h,i,j,k,l,m
      complex*16 echi(4),cechi(4)
      complex*16 epsi1(3,4),epsi(4),cepsi(4)
      double precision q5(4),q6(4),q7(4),q8(4),q9(4)

      include 'polvecs.f'
      include 'mom.f'
      include 'quarkonia.f'

      do i=1,3
         do j=1,3
            rho1chi(i,j)=wt(i)*wt(j+3)
         enddo
      enddo

cccccccccccccccc

      do i=1,4
         q7(i)=q(i,7)
         q6(i)=q(i,6)
         q5(i)=q(i,5)
      enddo

      call genpol1(7,epsi1)
      mvec=dsqrt(q(4,7)**2-q(3,7)**2-q(2,7)**2-q(1,7)**2)

      pnorm=2d0*sdot(q6,q7)**2*(mchi**2+mvec**2)/(mchi**2*mvec**2)/3d0

ccccccccccccccccccccc 
      
      do k=1,3
         do j=1,3
            
            do h=1,4
               epsi(h)=epsi1(k,h)
               cepsi(h)=conjg(epsi1(j,h))
            enddo
            
            call cdot(epsi,q6,ze7q6)
            call cdot(cepsi,q6,zce7q6)
            
            rho1psi(k,j)=0d0
            
            do l=1,3
               do h=1,3
                  
                  do m=1,4
                     echi(m)=echi1(l,m)            
                     cechi(m)=conjg(echi1(h,m))
                  enddo
   
                  call ccdot(echi,epsi,ze5e7)
                  call ccdot(cechi,cepsi,zce5ce7)
                  call ccdot(epsi,cepsi,zce7e7)
                  call ccdot(echi,cechi,zce5e5)
                  call cdot(echi,q6,ze5q6)
                  call cdot(cechi,q6,zce5q6)

                  rho1psi(k,j)=rho1psi(k,j)+(-zce5q6*ze5q6*zce7e7-
     &                 zce5e5*zce7q6*ze7q6+zce5q6*ze5e7
     &                 *zce7q6+ze5q6*zce5ce7*ze7q6)*rho1chi(l,h)

               enddo
            enddo     
         enddo
      enddo

ccccccccccc

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

      wtt=wt1/pnorm

      return
      end
