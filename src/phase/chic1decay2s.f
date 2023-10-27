ccc   spin correlations for chi_1 decay to 2 scalars
      subroutine chic1decay2s(wt,wtt)
      implicit none
      double precision wtt,wt1
      complex*16 zq7e5,zq7ce5,zq6ce5,zq6e5
      complex*16 wt(10)
      complex*16 rho1chi(3,3)
      integer h,i,j,k,l,m
      complex*16 echi(4),cechi(4)
      double precision q6(4),q7(4)

      include 'quarkonia.f'
      include 'polvecs.f'
      include 'mom.f'

      do i=1,3
         do j=1,3
            rho1chi(i,j)=wt(i)*wt(j+3)
         enddo
      enddo

cccccccccccccccc

      do k=1,4
         q6(k)=q(k,6)
         q7(k)=q(k,7)
      enddo

         wt1=0d0

        do l=1,3
           do h=1,3

              do m=1,4
                 echi(m)=echi1(l,m)
                 cechi(m)=conjg(echi1(h,m))
              enddo
              
              call cdot(echi,q6,zq6e5)
              call cdot(echi,q7,zq7e5)
              call cdot(cechi,q6,zq6ce5)
              call cdot(cechi,q7,zq7ce5)
            
              wt1=wt1+dble((zq6e5-zq7e5)*(zq6ce5-zq7ce5)
     &             /(mchi**2-4d0*m2b**2)*3d0*rho1chi(l,h))
 
           enddo
        enddo

        wtt=wt1

      return
      end
