ccc   integrates bare + screened amplitude over k_t
ccc   (two-photon induced processes)
      subroutine schimcgamion(p1x,p1y,p2x,p2y,outg)
      implicit none
      double precision x00p,wt,x00p2
      double precision xggmin,yp,ypmax,ypmin
      complex*16 zout,zout1,zoutg,zoutoff2
      double precision tpx,tpy,tp2,t12,t22,t11,phiq
      double precision sc,qtmax,qt,screeningionint,qtmin
      double precision p1xp,p2xp,p1yp,p2yp
      double precision hy,hqt,hphi,del,dbl
      double precision p1x,p1y,p2x,p2y
      integer p,i,nphi,jqt,jphi
      complex*16 out(4,10),x0(10),outg(10)
      complex*16 zouttest,zouttest1
      integer nk

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

      nphi=s2int*4
      nk=s2int*4

      del=0.5d0

      qtmax=0.5d0
      qtmin=0d0

      if(ionbreakup)then
         if(fAA.eq.'01')then
            xggmin=1d-4         ! 01
         else
            xggmin=1d-2         ! otherwise
         endif
      else
         xggmin=1d-2            ! otherwise
      endif

      ypmax=dlog(xggmin**2+qtmax**2)
      ypmin=dlog(xggmin**2+qtmin**2)

      hphi=2d0*pi/dble(nphi)

      hy=(ypmax-ypmin)/dble(nk)
      hqt=qtmax/dble(nk)

      zoutg=0d0

      do p=1,pol
         outg(p)=0d0
         do i=1,4
            out(i,p)=0d0
         enddo
      enddo

      call wtgengam

      if(sfac)then

         do jqt=1,nk

            yp=((ypmax-ypmin)*xikt4(jqt)+ypmax+ypmin)/2d0

            qt=dexp(yp)-xggmin**2
            tp2=qt
            qt=dsqrt(dabs(qt))

            sc=screeningionint(qt)


            zouttest=0d0
            zouttest1=0d0

            do jphi=1,nphi

               phiq=pi*(xiphi4(jphi)+1d0)

               tpx=qt*dcos(phiq)
               tpy=qt*dsin(phiq)

               wt=(tp2+xggmin**2)/2d0*(ypmax-ypmin)/2d0*pi
               wt=wt*wiphi4(jphi)*wikt4(jqt)

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

            call formfacgamoff_ion_surv(p,p1x,p1y,p2x,p2y,zoutoff2)

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

         else

            zout=-0.5d0*(ppa(p)+mma(p))*(p1x*p2x+p1y*p2y)
     &           -0.5d0*zi*(ppa(p)-mma(p))*(p1x*p2y-p2x*p1y)
     &           +0.5d0*(p1x*p2x-p1y*p2y+zi*(p1x*p2y+p1y*p2x))*mpa(p)
     &           +0.5d0*(p1x*p2x-p1y*p2y-zi*(p1x*p2y+p1y*p2x))*pma(p)

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

      return
      end
