ccc   spin correlations for chi_2 decay to 2 fermions
      subroutine chic2decay2f(wt,wtt)
      implicit none
      double precision q6q7,wt2,sdot,wtt,pnorm
      complex*16 wt(10)
      complex*16 rho2chi(5,5)
      integer h,i,j,k,l,m,o,p,n
      complex*16 echi(4,4),cechi(4,4)
      double precision q6(4),q7(4)
      double precision q6c(4),q7c(4)
      double precision g(4,4)

      include 'quarkonia.f'
      include 'polvecs.f'
      include 'mom.f'
      include 'vars.f'
      include 'partonmom4.f'

      do i=1,5
         do j=1,5
            rho2chi(i,j)=wt(i)*wt(j+5)
         enddo
      enddo

cccccccccccccccccccc

      wt2=0d0

      do k=1,4
         q6(k)=q(k,6)
         q7(k)=q(k,7)
      enddo

      q6q7=sdot(q6,q7)

      q6c(4)=q(4,6)
      q7c(4)=q(4,7)

      do j=1,3
         q6c(j)=-q(j,6)
         q7c(j)=-q(j,7)
      enddo

      q6c(4)=paa(4)
      q7c(4)=paa(4)
      do i=1,3
         q6c(i)=-paa(i)
         q7c(i)=paa(i)
      enddo

      do j=1,4
         do k=1,4
            g(j,k)=0d0
         enddo
      enddo

      g(1,1)=-1d0
      g(2,2)=-1d0
      g(3,3)=-1d0
      g(4,4)=1d0

      do m=1,5
         do n=1,5

            do h=1,4
               do l=1,4
                  echi(h,l)=echi2(m,h,l)
                  cechi(h,l)=conjg(echi2(n,h,l))
               enddo
            enddo

            do j=1,4
               do o=1,4
                  do k=1,4
                     do p=1,4
                        wt2=wt2+dble(4d0*echi(p,k)*cechi(o,j)*
     &                       (q6c(p)*q7c(o)+q6c(o)*q7c(p)
     &                       -mchi**2*g(p,o)/2d0)
     &                       *(q6c(k)-q7c(k))*(q6c(j)-q7c(j))
     &                       *rho2chi(m,n))
                     enddo
                  enddo
               enddo
            enddo

         enddo
      enddo

      pnorm=16d0*(mchi**2/4d0-m2b**2)*(mchi**2/2d0+4d0*m2b**2/3d0)/5d0
      wtt=wt2/pnorm

      return
      end
