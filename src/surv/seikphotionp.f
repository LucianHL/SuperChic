ccc   integrates bare + screened amplitude over k_t
ccc   (photoproduction processes)
      subroutine schimcphotionp(p1x,p1y,p2x,p2y,out)
      implicit none
      double precision del,x00p,x00pp,wt1,wt2,wt3,wt
      double precision tpx,tpy,tp2,t11,t12,t22
      double precision qtmax,qt,phiq
      double precision p1xp,p2xp,p1yp,p2yp,sc,screeningionint
      double precision p1x,p1y,p2x,p2y
      integer p,nphi,nqt,jqt,jphi
      complex*16 out(10),x0(10),out1(10),out2(10),x01(10)

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
      include 'gaussvars.f'

      del=0.5d0

      do p=1,pol
         out(p)=0d0
         out1(p)=0d0
         out2(p)=0d0
      enddo

      nphi=s2int*2
      nqt=s2int*2

      qtmax=0.5d0

      if(sfac)then

         do jqt=1,nqt

            qt=qtmax/2d0*(xikt(jqt)+1d0)

            tp2=qt**2

            sc=screeningionint(qt)

            do jphi=1,nphi

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

           do p=1,pol

              call formfacphotionp(1,t12,t22,x00p)
              call formfacphotionp(2,t12,t22,x00pp)

              x0(p)=x00p
              x01(p)=x00pp

              out(p)=out(p)+x0(p)*wt*sc*p2xp
              out1(p)=out1(p)+x0(p)*wt*sc*p2yp
              out2(p)=out2(p)+x01(p)*wt*sc

           enddo

      enddo
      enddo

      endif

      t11=p1x**2+p1y**2
      t22=p2x**2+p2y**2

      call formfacphotionp(1,t11,t22,x00p)
      call formfacphotionp(2,t11,t22,x00pp)

      do p=1,pol

         wt1=x00p*p2x
         wt2=x00p*p2y
         wt3=x00pp

         out(p)=out(p)+wt1
         out1(p)=out1(p)+wt2
c         out2(p)=out2(p)+wt3
         out2(p)=0d0   ! Remove as not physical form

         if(dabs(wt1).lt.cdabs(out(p))*del)out(p)=wt1
         if(dabs(wt2).lt.cdabs(out1(p))*del)out1(p)=wt2
c         if(dabs(wt3).lt.cdabs(out2(p))*del)out2(p)=wt3

      enddo

      do p=1,pol
         out(p)=dsqrt(cdabs(out(p))**2+cdabs(out1(p))**2
     &        +cdabs(out2(p))**2)
      enddo

      return
      end
