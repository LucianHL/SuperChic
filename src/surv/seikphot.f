ccc   integrates bare + screened amplitude over k_t
ccc   (photoproduction processes)
      subroutine schimcphot(p1x,p1y,p2x,p2y,out)
      implicit none
      double precision x00p,x00pp,wt
      double precision tpx,tpy,tp2,t22,t12,t11,sc,sc1
      double precision qt,phiq,p1xp,p1yp,p2xp,p2yp,hqt,hphi
      double precision p1x,p1y,p2x,p2y
      integer i1,i2,p,nphi,nqt,jqt,jphi
      complex*16 out(10),x0(10),out1(10),out2(10),x01(10)
      complex*16 screen(2,2)

      include 'nchan.f'
      include 'surv.f'
      include 'vars.f'
      include 'survpars.f'
      include 'polarization.f'
      include 'photo.f'
      include 'bpsi.f'
      include 'pi.f'
      include 'nsurv.f'
      include 'beam.f'

      do p=1,pol
         out(p)=0d0
         out1(p)=0d0
         out2(p)=0d0
      enddo

      nphi=s2int
      nqt=s2int*4

      hphi=2d0*pi/dble(nphi)
      hqt=2d0/dble(nqt)

      if(sfac)then

         do jqt=1,nqt

            qt=(dble(jqt)-0.5d0)*hqt

            tp2=qt**2

            do i1=1,nch
               do i2=1,nch
                  call screeningint(i1,i2,tp2,sc,sc1)
                  screen(i1,i2)=sc
               enddo
            enddo

            do jphi=1,nphi

               phiq=(dble(jphi)-0.5d0)*hphi

               tpx=qt*dcos(phiq)
               tpy=qt*dsin(phiq)
               wt=hphi*qt*hqt

               p1xp=p1x-tpx
               p1yp=p1y-tpy
               t12=p1xp**2+p1yp**2
               p2xp=tpx+p2x
               p2yp=tpy+p2y
               t22=p2xp**2+p2yp**2

           do p=1,pol

              do i1=1,nch
                 do i2=1,nch

                    call formfacphot(1,t12,t22,x00p)
                    call formfacphot(2,t12,t22,x00pp)

                    x0(p)=x00p*pp0(i1)*pp0(i2)/dble(nch)**2/gaa(i1)
     &                   /gaa(i2)
                    x01(p)=x00pp*pp0(i1)*pp0(i2)/dble(nch)**2/gaa(i1)
     &                   /gaa(i2)

                    if(prot.eq.1)then
                       out(p)=out(p)+x0(p)*wt*screen(i1,i2)*p1xp
                       out1(p)=out1(p)+x0(p)*wt*screen(i1,i2)*p1yp
                    else
                       out(p)=out(p)+x0(p)*wt*screen(i1,i2)*p2xp
                       out1(p)=out1(p)+x0(p)*wt*screen(i1,i2)*p2yp
                    endif

                    out2(p)=out2(p)+x01(p)*wt*screen(i1,i2)

                 enddo
              enddo

           enddo

      enddo
      enddo

      endif

      t11=p1x**2+p1y**2
      t22=p2x**2+p2y**2

      call formfacphot(1,t11,t22,x00p)
      call formfacphot(2,t11,t22,x00pp)

      if(prot.eq.1)then
         do p=1,pol
            out(p)=out(p)+x00p*p1x
            out1(p)=out1(p)+x00p*p1y
         enddo
      else
         do p=1,pol
            out(p)=out(p)+x00p*p2x
            out1(p)=out1(p)+x00p*p2y
         enddo
      endif

      do p=1,pol
         out2(p)=out2(p)+x00pp
      enddo

      do p=1,pol
         out(p)=dsqrt(cdabs(out(p))**2+cdabs(out1(p))**2
     &        +cdabs(out2(p))**2)
      enddo

      return
      end
