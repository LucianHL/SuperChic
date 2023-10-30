ccc   spin correlations for chi_2 decay to 2 scalars
      subroutine chic2decay2s(wt,wtt)
      implicit none
      double precision wt2,wt2a,wtt
      complex*16 wt(10)
      complex*16 rho2chi(5,5)
      complex*16 zwt2a,zwt2b
      integer h,i,j,k,l,m,o,p,n
      complex*16 echi(4,4),cechi(4,4)
      double precision q6c(4),q7c(4)

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
      
      
      wt2=0d0
      wt2a=0d0
  
   
              
      do m=1,5
         do n=1,5

            do h=1,4
               do l=1,4
                  echi(h,l)=echi2(m,h,l)
                  cechi(h,l)=conjg(echi2(n,h,l))
               enddo
            enddo

            zwt2a=(0d0,0d0)
            zwt2b=(0d0,0d0)
            
            do j=1,4
               do k=1,4
                  do p=1,4
                     do o=1,4
                    wt2=wt2+dble(echi(j,k)*cechi(p,o)*(q6c(j)-q7c(j))
     &                       *(q6c(k)-q7c(k))*(q6c(p)-q7c(p))
     &                       *(q6c(o)-q7c(o))*rho2chi(m,n))
                     enddo
                  enddo
               enddo
            enddo

            do j=1,4
               do k=1,4
                  zwt2a=zwt2a+echi(j,k)*(q6c(j)-q7c(j))
     &                 *(q6c(k)-q7c(k))
                  zwt2b=zwt2b+cechi(j,k)*(q6c(j)-q7c(j))
     &                 *(q6c(k)-q7c(k))
               enddo
            enddo
                  
            wt2a=wt2a+dble(zwt2a*zwt2b*rho2chi(m,n))

            
            
         enddo
      enddo

      wtt=wt2/(mchi**2-4d0*m2b**2)**2*15d0/2d0

      
      return
      end
