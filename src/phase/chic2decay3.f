ccc   spin correlations for chi_2 decay to mu+mu- + gamma via intermediate
ccc   vector meson
      subroutine chic2decay3(wt,wtt)
      implicit none
      double precision q8q9,q6q7,sdot,wt1,pnorm,wtt,mvec,mxx
      complex*16 zq6epsi1p,ze7q8,zce7q8,zq6cepsi1,zce7e7
      complex*16 wt(10)
      complex*16 rho1psi(3,3),rho2chi(5,5)
      integer h,i,j,k,l,m,o,p,r,d,n
      complex*16 epsi1(3,4),epsi(4),cepsi(4)
      complex*16 epsic(4),cepsic(4)
      double precision q6(4),q7(4),q8(4),q9(4)
      double precision q6c(4),q7c(4)
      double precision g(4,4)
      double precision pv(4)
      complex*16 echi(4,4),cechi(4,4)

      include 'quarkonia.f'
      include 'polvecs.f'
      include 'mom.f'
      include 'partonmom4.f'
      
      do i=1,5
         do j=1,5
            rho2chi(i,j)=wt(i)*wt(j+5)
         enddo
      enddo

ccccccccccccccccccc

      mvec=dsqrt(q(4,7)**2-q(3,7)**2-q(2,7)**2-q(1,7)**2)
      mxx=dsqrt(q(4,5)**2-q(3,5)**2-q(2,5)**2-q(1,5)**2)
      
      pv(4)=mxx-paa(4)
      do i=1,3
         pv(i)=-paa(i)
      enddo

      call genpol1rf(pv,epsi1)
      
      do j=1,4
         do k=1,4
            g(j,k)=0d0
         enddo
      enddo

      g(1,1)=-1d0
      g(2,2)=-1d0
      g(3,3)=-1d0
      g(4,4)=1d0

cccccccccccccccccccc

c     generate contravariant 4-vectors

       q6c(4)=paa(4)
       q7c(4)=pv(4)
       do j=1,3
          q6c(j)=-paa(j)
          q7c(j)=-pv(j)
       enddo

       do j=1,4
          q6(j)=paa(j)
          q7(j)=pv(j)
          q8(j)=pbb(j)
          q9(j)=-pbb(j)
       enddo
       q9(4)=pbb(4)
       
       q6q7=sdot(q6,q7)

ccccccccccccccccc

c     normalisation

      do k=1,3
         do j=1,3

            do h=1,4
               cepsi(h)=conjg(epsi1(k,h))
               epsi(h)=epsi1(j,h)
            enddo
            
            epsic(4)=epsi(4)
            cepsic(4)=cepsi(4)
            do m=1,3
               epsic(m)=-epsi(m)
               cepsic(m)=-cepsi(m)
            enddo
            
            call cdot(cepsi,q6,zq6cepsi1)
            call cdot(epsi,q6,zq6epsi1p)

            rho1psi(k,j)=0d0
            
            do m=1,5
               do n=1,5

                  do h=1,4
                     do l=1,4
                        echi(h,l)=echi2(m,h,l)
                        cechi(h,l)=conjg(echi2(n,h,l))
                     enddo
                  enddo
                  
                  do o=1,4
                     do p=1,4
                        do d=1,4
                           do r=1,4
                              
          rho1psi(k,j)=rho1psi(k,j)+(echi(r,d)*
     &    cechi(p,o)*q6c(r)*(q6q7*epsic(o)-zq6epsi1p*q7c(o)
     &    )*(cepsic(d)*q7c(p)-cepsic(p)*q7c(d))+
     &    cechi(r,d)*echi(p,o)*q6c(r)*(q6q7*cepsic(o)-zq6cepsi1
     &    *q7c(o))*(epsic(d)*q7c(p)-epsic(p)*q7c(d))
     &    -echi(r,d)*cechi(p,o)*(q6c(p)*q6c(r)*(mvec**2*epsic(o)
     &    *cepsic(d)-q7c(d)*q7c(o))+g(r,p)*(q6q7*epsic(o)-
     &    zq6epsi1p*q7c(o))*(q6q7*cepsic(d)
     &    -zq6cepsi1*q7c(d))))*rho2chi(m,n)
            
                           enddo
                        enddo
                     enddo
                  enddo
      
               enddo
            enddo
            
         enddo
      enddo
         
      pnorm=10d0*(1d0-q6q7/mchi**2)*q6q7**2/15d0

ccccccccccccccccccccccccccc

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

            wt1=wt1+dble((-(q8q9+mmu**2)*zce7e7-2d0*ze7q8*
     &zce7q8)*rho1psi(k,j)/(2d0*(q8q9+2d0*mmu**2))*3d0)

         enddo
      enddo

      wtt=wt1/pnorm

      return
      end
