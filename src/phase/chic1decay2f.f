ccc   spin correlations for chi_1 decay to 2 fermions
      subroutine chic1decay2f(wt,wtt)
      implicit none
      double precision q6q7,wt1,sdot,wtt
      complex*16 zq6ce5,zq6e5,ze5ce5
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
              call cdot(cechi,q6,zq6ce5)
              call ccdot(echi,cechi,ze5ce5)
              q6q7=sdot(q6,q7)
              
         wt1=wt1+(-(q6q7+m2b**2)*ze5ce5-2d0*zq6e5*
     &      zq6ce5)*rho1chi(l,h)/(2d0*(q6q7+2d0*m2b**2))*3d0   

           enddo
        enddo

        wtt=wt1

      return
      end
