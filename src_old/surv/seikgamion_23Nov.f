ccc   integrates bare + screened amplitude over k_t
ccc   (two-photon induced processes)
      subroutine schimcgamion(p1x,p1y,p2x,p2y,outg)
      implicit none
      double precision x00p,wt,x00p2
      double precision xggmin,yp,ypmax,ypmin,yqmin,yqmax,yq
      complex*16 zout,zout1,zoutg,zoutoff2
      double precision tpx,tpy,tp2,t12,t22,t11,phiq
      double precision sc,qtmax,qt,screeningionint,qtmin
      double precision p1xp,p2xp,p1yp,p2yp
      double precision hy,hqt,hphi,del,dbl
      double precision hy2,y2max,y2min,y2
      double precision p1x,p1y,p2x,p2y
      integer jx,jy,i1,i2,p,i,nphi,nqt,jqt,jphi
      complex*16 out(4,10),x0(10),x00(10),outg(10)
      complex*16 screen(2,2),zouttest,zouttest1
      integer nni,nk
      complex*16 outt(0:200,0:10)

      include 'ppamp.f'
      include 'nchan.f'
      include 'surv.f'
      include 'vars.f'
      include 'survpars.f'
      include 'polarization.f'
      include 'photo.f'
      include 'bpsi.f'
      include 'proc.f'
      include 'zi.f'
      include 'mandelstam.f'
      include 'pi.f'
      include 'nsurv.f'
      include 'inparticle.f'
      include 'zarr.f'
      include 'beam.f'
      include 'gaussvars.f'
      include 'p0Xn.f'
      include 'diss.f'
      include 'eff.f'

      if(neff.eq.1)then
         write(404,*)'neffs2 = 1'
         call flush(404)
      endif

      nphi=s2int*4
      nk=s2int*4

c      nphi=100
c      nk=100

      del=0.5d0

      qtmax=0.5d0
c      qtmax=1d0


c      qtmax=0.5d0
c      qtmax=2d-3
      qtmin=1d-3
      qtmin=0d0

      xggmin=1d-2

c      xggmin=0d0

      ypmax=dlog(xggmin**2+qtmax**2)
      ypmin=dlog(xggmin**2+qtmin**2)

      hphi=2d0*pi/dble(nphi)

      y2max=-1d0/(qtmax**2+xggmin**2)
      y2min=-1d0/(qtmin**2+xggmin**2)
      hy2=(y2max-y2min)/dble(nk)

      y2max=-1d0/(qtmax+xggmin)
      y2min=-1d0/(qtmin+xggmin)


c      nphi=s2int*16
c      nk=s2int*16


      yqmax=dlog(xggmin+qtmax)
      yqmin=dlog(xggmin+qtmin)

      hy=(ypmax-ypmin)/dble(nk)
      hqt=qtmax/dble(nk)

      zoutg=0d0

c      do i1=0,nk
c         do p=0,10
c            outt(i1,p)=0d0
c         enddo
c      enddo

      do p=1,pol
         outg(p)=0d0
         do i=1,4
            out(i,p)=0d0
         enddo
      enddo

      call wtgengam

      if(sfac)then

         do jqt=1,nk

cc            yp=((ypmax-ypmin)*xikt(jqt)+ypmax+ypmin)/2d0
            yp=ypmin+(ypmax-ypmin)*(dble(jqt)-0.5d0)/dble(nk)

            qt=dexp(yp)-xggmin**2
            tp2=qt
            qt=dsqrt(dabs(qt))

c            y2=((y2max-y2min)*xikt(jqt)+y2max+y2min)/2d0
            y2=y2min+(y2max-y2min)*(dble(jqt)-0.5d0)/dble(nk)
c            qt=-1d0/y2-xggmin**2
c            qt=dsqrt(qt)

c            yq=yqmin+(yqmax-yqmin)*(dble(jqt)-0.5d0)/dble(nk)
c             yq=((yqmax-yqmin)*xikt(jqt)+yqmax+yqmin)/2d0
            y2=((y2max-y2min)*xikt(jqt)+y2max+y2min)/2d0
c            qt=((qtmax-qtmin)*xikt(jqt)+qtmax+qtmin)/2d0
c            qt=dsqrt(qt)
c            qt=dexp(yq)-xggmin

            qt=-1d0/y2-xggmin

c            qt=qtmin+(qtmax-qtmin)*(dble(jqt)-0.5d0)/dble(nk)

            sc=screeningionint(qt)

            zouttest=0d0
            zouttest1=0d0

            do jphi=1,nphi

               phiq=pi*(xiphi(jphi)+1d0)
c               phiq=2d0*pi*dble(jphi-1)/dble(nphi)

               tpx=qt*dcos(phiq)
               tpy=qt*dsin(phiq)

c               wt=(tp2+xggmin**2)*(ypmax-ypmin)*2d0*pi
c     &              /dble(nphi)/dble(nk)/2d0

c              wt=(qt**2+xggmin**2)**2*(y2max-y2min)*2d0*pi
c     &              /dble(nphi)/dble(nk)/2d0

c              wt=(qt+xggmin)**2*(y2max-y2min)*2d0*pi
c     &              /dble(nphi)/dble(nk)*qt

c              wt=(qt+xggmin)*(yqmax-yqmin)*2d0*pi
c     &              /dble(nphi)/dble(nk)*qt

cc               wt=(tp2+xggmin**2)/2d0*(ypmax-ypmin)/2d0*pi
cc               wt=1d0/2d0*(qtmax-qtmin)/2d0*pi*qt*2d0
                wt=(qt+xggmin)**2/2d0*(y2max-y2min)/2d0*pi*qt*2d0
cc               wt=(qt**2+xggmin**2)**2/2d0*(y2max-y2min)/2d0*pi

               wt=wt*wiphi(jphi)*wikt(jqt)

cc               wt=(qtmax-qtmin)/dble(nk)*2d0*qt*2d0*pi/dble(nphi)

               p1xp=p1x-tpx
               p1yp=p1y-tpy
               t12=p1xp**2+p1yp**2
               p2xp=tpx+p2x
               p2yp=tpy+p2y
               t22=p2xp**2+p2yp**2

               if(beam.eq.'ion')then
                  if(offshell)then
                  else
                     call formfacgamion(t12,t22,x00p)
                  endif
               elseif(beam.eq.'ionp')then
                  call formfacgamionp(1,t12,t22,x00p)
                  call formfacgamionp(2,t12,t22,x00p2)
               endif

           do p=1,pol

              if(offshell)then

                 call formfacgamoff_ion_surv(p,p1xp,p1yp,p2xp,p2yp
     &                ,zout1)
                 zout=zout1
                 outg(p)=outg(p)+wt*sc*zout


c                 if(qt.gt.0.1d0)wt=0d0

                 zouttest=zouttest+wt*sc*zout
                 zouttest1=zouttest1+wt*sc

              else

                 x0(p)=x00p

                 zout=-0.5d0*(ppa(p)+mma(p))*(p1xp*p2xp+p1yp*p2yp)
     &                -0.5d0*zi*(ppa(p)-mma(p))*(p1xp*p2yp-p2xp*p1yp)
     &                +0.5d0*(p1xp*p2xp-p1yp*p2yp
     &                +zi*(p1xp*p2yp+p1yp*p2xp))*mpa(p)
     &                +0.5d0*(p1xp*p2xp-p1yp*p2yp
     &                -zi*(p1xp*p2yp+p1yp*p2xp))*pma(p)

                 zout=zout*2d0

                 outg(p)=outg(p)+x0(p)*wt
     &                *sc*zout

                 zouttest=zouttest+wt*sc*zout*x0(p)
                 zouttest1=zouttest1+wt*sc*x0(p)


                 if(p.eq.1)then
                    zoutg=zoutg+x00p2*wt*sc
                 endif

              endif


           enddo

        enddo

c        print*,yp,qt,cdabs(zouttest),cdabs(zouttest1)

      enddo

      endif

      t11=p1x**2+p1y**2
      t22=p2x**2+p2y**2

      if(beam.eq.'ion')then
         call formfacgamion(t11,t22,x00p)
         x00p2=0d0
      elseif(beam.eq.'ionp')then
         call formfacgamionp(1,t11,t22,x00p)
         call formfacgamionp(2,t11,t22,x00p2)
      endif

      do p=1,pol

         if(offshell)then

c            call formfacgamoff_ion(p,p1x,p1y,p2x,p2y,zoutoff2s)
            call formfacgamoff_ion_surv(p,p1x,p1y,p2x,p2y,zoutoff2)

c            print*,cdabs(zoutoff2s)/cdabs(zoutoff2)

c            if(sfac)then
c               call formfacgamoff_ion_surv(p,p1x,p1y,p2x,p2y,zoutoff2s)
c            else
c               outg(p)=dsqrt(cdabs(zoutoff2)**2)
c            endif



c            outg(p)=outg(p)+zoutoff2


            if(ionbreakup)then
               if(sfac)then
                  if(fAA.eq.'00')then ! exp(-omega)*p and not 1-exp(-omega)*p (more stable) for other cases
                     outg(p)=outg(p)+zoutoff2
                  endif
               else
                  outg(p)=outg(p)+zoutoff2
               endif
            else
               outg(p)=outg(p)+zoutoff2
            endif

            outg(p)=dsqrt(cdabs(outg(p))**2)


c            outg(p)=outg(p)*cdabs(zoutoff2s)/cdabs(zoutoff2)

c            print*,p,outg(p)

         else

            zout=-0.5d0*(ppa(p)+mma(p))*(p1x*p2x+p1y*p2y)
     &           -0.5d0*zi*(ppa(p)-mma(p))*(p1x*p2y-p2x*p1y)
     &           +0.5d0*(p1x*p2x-p1y*p2y+zi*(p1x*p2y+p1y*p2x))*mpa(p)
     &           +0.5d0*(p1x*p2x-p1y*p2y-zi*(p1x*p2y+p1y*p2x))*pma(p)

c            zout=(p1x*p2x+p1y*p2y)


            if(ionbreakup)then
               if(sfac)then
                  if(fAA.eq.'00')then ! exp(-omega)*p and not 1-exp(-omega)*p (more stable) for other cases
                     outg(p)=outg(p)+zout*x00p*2d0
                  endif
               else
                  outg(p)=outg(p)+zout*x00p*2d0
               endif
            else
               outg(p)=outg(p)+zout*x00p*2d0
            endif




            if(cdabs(zout*x00p*2d0).lt.cdabs(outg(p))*del)
     &           outg(p)=zout*x00p*2d0

            dbl=cdabs(outg(p))
            if(dabs(x00p2).lt.cdabs(zoutg)*del)zoutg=x00p2
            dbl=dsqrt(dbl**2+cdabs(zoutg+x00p2)**2*pincarr(p))
            outg(p)=dbl

         endif

      enddo

      if(neff.eq.1)then
         write(404,*)'neffs2 = 1 a'
         call flush(404)
      endif

      return
      end
