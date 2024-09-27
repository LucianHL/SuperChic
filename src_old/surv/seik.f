ccc   integrates bare + screened amplitude over k_t
ccc   (QCD induced processes)
      subroutine schimc(p1x,p1y,p2x,p2y,out)
      implicit none
      double precision x00p,mu
      double precision p1x,p1y,p2x,p2y,tpx,tpy,tp2,t11,t22,t12
      double precision sc,sc1
      double precision qtmax,qt,phiq,hphi,hqt
      double precision p1xp,p2xp,p1yp,p2yp
      double precision wt
      integer i1,i2,p,nphi,nqt,jqt,jphi,nnphi
      complex*16 out(10),x0(10),x00(10),x0p(10)
      complex*16 screen(2,2)
      complex*16 outt(0:200,10)

      include 'nchan.f'
      include 'mt.f'
      include 'surv.f'
      include 'vars.f'
      include 'survpars.f'
      include 'polarization.f'
      include 'nsurv.f'
      include 'pi.f'
      include 'gaussvars.f'

      call setmu(mu)

      do p=1,pol
         out(p)=0d0
      enddo

      nphi=s2int*2
      nqt=s2int*2

      qtmax=1.5d0

      hphi=2d0*pi/dble(nphi)
      hqt=qtmax/dble(nqt)

      do i1=0,nqt
         do p=1,10
            outt(i1,p)=0d0
         enddo
      enddo

      if(sfac)then

        do jqt=1,nqt

           qt=qtmax/2d0*(xikt(jqt)+1d0)

            tp2=qt**2

            do i1=1,nch
               do i2=1,nch
                  call screeningint(i1,i2,tp2,sc,sc1)
                  screen(i1,i2)=sc
               enddo
            enddo

            nnphi=nphi

           do jphi=1,nnphi

              phiq=pi*(xiphi(jphi)+1d0)

               tpx=qt*dcos(phiq)
               tpy=qt*dsin(phiq)
               wt=qt*wiphi(jphi)*wikt(jqt)*qtmax/2d0*pi

            p1xp=p1x-tpx
            p1yp=p1y-tpy
            t12=p1xp**2+p1yp**2
            p2xp=tpx+p2x
            p2yp=tpy+p2y
            t22=p2xp**2+p2yp**2

          call bare(mu,p1xp,p1yp,p2xp,p2yp,x0p)

           do p=1,pol

              do i1=1,nch
                 do i2=1,nch

                    x0(p)=x0p(p)*pp0(i1)*pp0(i2)/dble(nch)**2
                    x0(p)=x0(p)*dexp(-((t12+0.08d0+bb0(i1))*bex(i1))
     &                   **cc0(i1)+(bex(i1)*(bb0(i1)+0.08d0))**cc0(i1))
                    x0(p)=x0(p)*dexp(-((t22+0.08d0+bb0(i2))*bex(i2))
     &                   **cc0(i2)+(bex(i2)*(bb0(i2)+0.08d0))**cc0(i2))

                    out(p)=out(p)+x0(p)*wt*screen(i1,i2)

                 enddo
              enddo

           enddo

      enddo
      enddo

      endif

      t11=p1x**2+p1y**2
      t22=p2x**2+p2y**2

      call formfac(t11,t22,x00p)
      call bare(mu,p1x,p1y,p2x,p2y,x00)


      do p=1,pol
         out(p)=out(p)+x00(p)*x00p
      enddo

      return
      end
