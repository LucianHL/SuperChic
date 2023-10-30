ccc   spin correlations for j_psi decay to mu+mu-
      subroutine jpsidecay(wt,wtt)
      implicit none
      double precision wtt,q10q11,q8q9,sdot,sh,pnorm
      complex*16 ze6q8,zce7e7,ze7q10,zce7q10,zce6q8,
     & zce6e6
      complex*16 wt(10),wt1
      complex*16 rho1psi(4,4)
      integer h,i,j,k
      complex*16 epsi1(3,4),epsia(4),cepsia(4)
      complex*16 epsi2(3,4),epsib(4),cepsib(4)
      double precision pcm(4),pboo(4),plb(4)
      double precision q8(4),q9(4),q10(4),q11(4)

      include 'polvecs.f'
      include 'mom.f'
      include 'quarkonia.f'
      include 'polarization.f'
      include 'partonmom2.f'
      include 'partonmom4.f'

      do i=1,pol
         do j=1,pol
            rho1psi(i,j)=wt(i)*conjg(wt(j))
         enddo
      enddo

cccccccccccccccccccc

      do k=1,4
         q8(k)=q(k,8)
         q9(k)=q(k,9)
         q10(k)=q(k,10)
         q11(k)=q(k,11)
      enddo

      q8q9=sdot(q8,q9)
      q10q11=sdot(q10,q11)

      do k=1,4
         q8(k)=paa(k)
         q10(k)=pbb(k)
      enddo

      do k=1,4
       q(k,16)=p1(k)
       q(k,17)=p2(k)
      enddo
      
      call genpol1(16,epsi1)
      call genpol1(17,epsi2)

      do k=1,3
         pboo(k)=-p1(k)
      enddo
      pboo(4)=p1(4)
      sh=dsqrt(pboo(4)**2-pboo(3)**2-pboo(2)**2-pboo(1)**2)

      do k=1,4
         pcm(k)=paa(k)
      enddo
      call boost(sh,pboo,pcm,plb) 
      do k=1,4
         q8(k)=plb(k)
      enddo

ccccccccccccc

      do k=1,3
         pboo(k)=-p2(k)
      enddo
      pboo(4)=p2(4)
      sh=dsqrt(pboo(4)**2-pboo(3)**2-pboo(2)**2-pboo(1)**2)

      do k=1,4
         pcm(k)=pbb(k)
      enddo
      call boost(sh,pboo,pcm,plb) 
      do k=1,4
         q10(k)=plb(k)
      enddo

ccccccc
   
      wt1=(0d0,0d0)
      pnorm=0d0
      
      do j=1,pol
         do k=1,pol

            do h=1,4
               if(j.eq.1)then
                  cepsia(h)=conjg(epsi1(1,h))
                  cepsib(h)=conjg(epsi2(2,h))
               elseif(j.eq.2)then
                  cepsia(h)=conjg(epsi1(1,h))
                  cepsib(h)=conjg(epsi2(1,h))
               elseif(j.eq.3)then
                  cepsia(h)=conjg(epsi1(3,h))
                  cepsib(h)=conjg(epsi2(3,h))
               elseif(j.eq.4)then
                  cepsia(h)=conjg(epsi1(1,h))
                  cepsib(h)=conjg(epsi2(3,h))
               endif
               if(k.eq.1)then
                  epsia(h)=epsi1(1,h)
                  epsib(h)=epsi2(2,h)
               elseif(k.eq.2)then
                  epsia(h)=epsi1(1,h)
                  epsib(h)=epsi2(1,h)
               elseif(k.eq.3)then
                  epsia(h)=epsi1(3,h)
                  epsib(h)=epsi2(3,h)
               elseif(k.eq.4)then
                  epsia(h)=epsi1(1,h)
                  epsib(h)=epsi2(3,h)
               endif
            enddo

               call cdot(epsia,q8,ze6q8)
               call cdot(cepsia,q8,zce6q8)
               call ccdot(epsia,cepsia,zce6e6)
               call cdot(epsib,q10,ze7q10)
               call cdot(cepsib,q10,zce7q10)
               call ccdot(epsib,cepsib,zce7e7)

            wt1=wt1+(-(q8q9+mmu**2)*zce6e6-2d0*ze6q8*
     &zce6q8)*rho1psi(k,j)/(2d0*(q8q9+2d0*mmu**2))*3d0
     &*(-(q10q11+mmu**2)*zce7e7-2d0*ze7q10*
     &zce7q10)/(2d0*(q10q11+2d0*mmu**2))*3d0

         enddo
      enddo

      wtt=dble(wt1)

      return
      end
