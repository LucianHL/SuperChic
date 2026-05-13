ccc   integrates bare + screened amplitude over k_t
ccc   (two-photon induced processes)
      subroutine schimcgamion(p1x,p1y,p2x,p2y,outg)
      implicit none
      double precision x00p,wt,x00p2
      double precision xggmin,yp,ypmax,ypmin
      complex*16 zout,zout1,zoutg,zoutoff2,zoutoff2_nos2
      double precision tpx,tpy,tp2,t12,t22,t11,phiq
      double precision sc,qtmax,qt,screeningionint,qtmin
      double precision p1xp,p2xp,p1yp,p2yp
      double precision hy,hqt,hphi,del,dbl
      double precision p1x,p1y,p2x,p2y,denom,wtt_nos2,wtt_s2
      integer p,i,nphi,jqt,jphi
      complex*16 out(4,10),x0(10),outg(10),outg2(10)
      complex*16 zouttest,zouttest1,zout2,zoutoff2_2
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
      include 'ion_inel.f'

      nphi=s2int*4
      nk=s2int*4

      nphi=24*4

      wtt_s2=0d0
      wtt_nos2=0d0

      del=0.1d0

      qtmin=1d-7   
      qtmax=0.5d0

      if(ionbreakup)then
         if(int_01)then
            xggmin=1d-5         ! 01
            qtmax=0.1d0
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
         outg2(p)=0d0
         do i=1,4
            out(i,p)=0d0
         enddo
      enddo

      call wtgengam

      if(sfac)then

         do jqt=1,nk

           yp=((ypmax-ypmin)*xikt4(jqt)+ypmax+ypmin)/2d0
c            yp=ypmin+(ypmax-ypmin)*(dble(jqt)-0.5d0)/dble(nk)

c            yp=ypmin

            qt=dexp(yp)-xggmin**2
c            qt=xggmin*dtan(xggmin*yp/2d0)
c            qt=qt**2

c            qt=(-2d0/3d0/yp)**(1d0/3d0)
c            qt=qt**2

            
            tp2=qt
            qt=dsqrt(dabs(qt))

c            print*,qt

            sc=screeningionint(qt)

c            print*,qt,sc


            zouttest=0d0
            zouttest1=0d0

            do jphi=1,nphi

               phiq=pi*(xiphi4(jphi)+1d0)
c               phiq=(dble(jphi)-0.5d0)*2d0*pi/dble(nphi)

               tpx=qt*dcos(phiq)
               tpy=qt*dsin(phiq)

               wt=(tp2+xggmin**2)/2d0*(ypmax-ypmin)/2d0*pi
               wt=wt*wiphi4(jphi)*wikt4(jqt)

c               wt=(tp2+xggmin**2)/2d0
c               wt=wt*2d0*pi/dble(nphi)*(ypmax-ypmin)/dble(nk)

c               wt=wt*qt               
c               wt=qt**3/2d0


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
                  if(ion_inel)then
                  call formfacgamion_inel(1,t12,t22,x00p)
                  call formfacgamion_inel(2,t12,t22,x00p2)
                  else
                  call formfacgamionp(1,t12,t22,x00p)
                  call formfacgamionp(2,t12,t22,x00p2)
                  endif
               endif

           do p=1,pol  ! LHL fixpar

              if(offshell)then

               if(ion_inel)then
         call formfacgamoff_ionp_surv(1,p,p1xp,p1yp,p2xp,p2yp,zout1)
         call formfacgamoff_ionp_surv(2,p,p1xp,p1yp,p2xp,p2yp,zout2)
         outg2(p)=outg2(p)+wt*sc*zout2

               else
               call formfacgamoff_ion_surv(p,p1xp,p1yp,p2xp,p2yp,zout1)
               endif

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


                  if(jphi.eq.1)then
c                  print*,qt,sc,wt*sc*zout*x0(p)
                  endif

                 if(p.eq.1)then
                    zoutg=zoutg+x00p2*wt*sc
                 endif

              endif


           enddo

        enddo

      enddo

      endif

c      if(sfac)then
c      do p=1,pol
c      print*,outg(p)
c      enddo
c      endif

      t11=p1x**2+p1y**2
      t22=p2x**2+p2y**2

      if(beam.eq.'ion')then
            call formfacgamion(t11,t22,x00p)
           x00p2=0d0
      elseif(beam.eq.'ionp')then
         if(ion_inel)then
         if(ion_em)then
         call formfacgamion_ion_em(1,t11,t22,x00p)
         call formfacgamion_ion_em(2,t11,t22,x00p2)
         else
         call formfacgamion_inel(1,t11,t22,x00p)
         call formfacgamion_inel(2,t11,t22,x00p2)
         endif
         else
         call formfacgamionp(1,t11,t22,x00p)
         call formfacgamionp(2,t11,t22,x00p2)
         endif
      endif

      do p=1,pol  ! LHL fixpar

         if(ion_em)then


         call formfacgamoff(p,p1x,p1y,p2x,p2y,zoutoff2)
         outg(p)=dsqrt(cdabs(zoutoff2)**2)
      
         elseif(offshell)then

            if(ion_inel)then
            if(sfac)then
            call formfacgamoff_ionp_surv(1,p,p1x,p1y,p2x,p2y,zoutoff2)
            call formfacgamoff_ionp_surv(2,p,p1x,p1y,p2x,p2y,zoutoff2_2)
            call formfacgamoff_ionp(p,p1x,p1y,p2x,p2y,zoutoff2_nos2)
            else
            call formfacgamoff_ionp(p,p1x,p1y,p2x,p2y,zoutoff2)
            endif
            else
            if(sfac)then
            call formfacgamoff_ion_surv(p,p1x,p1y,p2x,p2y,zoutoff2)
            call formfacgamoff_ion(p,p1x,p1y,p2x,p2y,zoutoff2_nos2)
            else
            call formfacgamoff_ion(p,p1x,p1y,p2x,p2y,zoutoff2)
c            call formfacgamoff_ion_surv(p,p1x,p1y,p2x,p2y,zoutoff2)
            endif
            endif

            if(ionbreakup)then
               if(sfac)then
               if(wrho)then
               else
                  if(fAA.eq.'00'.or.fAA.eq.'AA'.or.fAA.eq.'A0')then ! exp(-omega)*p and not 1-exp(-omega)*p (more stable) for other cases
                     outg(p)=outg(p)+zoutoff2
                     outg2(p)=outg2(p)+zoutoff2_2
                  endif
               endif
               else
                  outg(p)=outg(p)+zoutoff2
                  outg2(p)=outg2(p)+zoutoff2_2
               endif
            else
               if(sfac)then
               endif
               outg(p)=outg(p)+zoutoff2
               outg2(p)=outg2(p)+zoutoff2_2
            endif

            outg(p)=dsqrt(cdabs(outg(p))**2+cdabs(outg2(p)))
            if(sfac)then

            denom=dsqrt(cdabs(zoutoff2)**2+cdabs(zoutoff2_2)**2)


            if(denom.eq.0d0)then
            outg(p)=0d0
            else
            outg(p)=outg(p)*cdabs(zoutoff2_nos2)/denom
c            outg(p)=cdabs(zoutoff2_nos2)
            endif
            endif

         else

            zout=-0.5d0*(ppa(p)+mma(p))*(p1x*p2x+p1y*p2y)
     &           -0.5d0*zi*(ppa(p)-mma(p))*(p1x*p2y-p2x*p1y)
     &           +0.5d0*(p1x*p2x-p1y*p2y+zi*(p1x*p2y+p1y*p2x))*mpa(p)
     &           +0.5d0*(p1x*p2x-p1y*p2y-zi*(p1x*p2y+p1y*p2x))*pma(p)


c            print*,p1x,p2x,p1y,p2y
c            print*,p1x*p2x-p1y*p2y+zi*(p1x*p2y+p1y*p2x)
c            print*,zout
c            print*,mpa(p),pma(p)
c            print*,''

c            if(sfac)print*,p,zout*x00p*2d0,outg(p),outg(p)+zout*x00p*2d0
c            if(sfac)wtt_nos2=wtt_nos2+cdabs(zout*x00p*2d0)**2
c            if(sfac)wtt_s2=wtt_s2+cdabs(outg(p)+zout*x00p*2d0)**2

            if(ionbreakup)then
               if(sfac)then
               if(wrho)then
               else
                  if(fAA.eq.'00'.or.fAA.eq.'AA'.or.fAA.eq.'A0')then ! exp(-omega)*p and not 1-exp(-omega)*p (more stable) for other cases
                     outg(p)=outg(p)+zout*x00p*2d0
                  endif
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

      subroutine schimcgamion_s2int(p1x,p1y,p2x,p2y,outg)
      implicit none
      double precision p1x,p1y,p2x,p2y,wtt,wtsum,wtt1,wtt2
      double precision del,diff,wtt3
      complex*16 outg(10)
      integer itvar
      common/itv/itvar

      include 'nsurv.f'
      include 'polarization.f'


      itvar=itvar+1

      del=1d-2

      call schimcgamion_s2var(16,p1x,p1y,p2x,p2y,outg)
      wtt1=wtsum(outg)

c      return

      if(dabs(wtt1).lt.1d-50)return

      call schimcgamion_s2var(32,p1x,p1y,p2x,p2y,outg)
      wtt2=wtsum(outg)

      diff=dabs(wtt1/wtt2-1d0)

      if(diff.lt.del)return

      call schimcgamion_s2var(64,p1x,p1y,p2x,p2y,outg)
      wtt3=wtsum(outg)
      diff=dabs(wtt2/wtt3-1d0)

      if(diff.lt.del)return

      wtt2=wtt3

      call schimcgamion_s2var(128,p1x,p1y,p2x,p2y,outg)
      wtt3=wtsum(outg)
      diff=dabs(wtt2/wtt3-1d0)

      if(diff.lt.del)return

      wtt2=wtt3

      call schimcgamion_s2var(256,p1x,p1y,p2x,p2y,outg)
      wtt3=wtsum(outg)
      diff=dabs(wtt2/wtt3-1d0)

      if(diff.lt.del)return

      wtt2=wtt3

      call schimcgamion_s2var(512,p1x,p1y,p2x,p2y,outg)
      wtt3=wtsum(outg)
      diff=dabs(wtt2/wtt3-1d0)

      if(diff.lt.del)return

      wtt2=wtt3

      call schimcgamion_s2var(256,p1x,p1y,p2x,p2y,outg)
      wtt3=wtsum(outg)
      diff=dabs(wtt2/wtt3-1d0)

      return
      end

      subroutine schimcgamion_s2int1(p1x,p1y,p2x,p2y,outg)
      implicit none
      double precision p1x,p1y,p2x,p2y,wtt,wtsum,wtt1,wtt2
      double precision del,diff,wtt3
      complex*16 outg(10)

      include 'nsurv.f'
      include 'polarization.f'

      call schimcgamion_s2var(64,p1x,p1y,p2x,p2y,outg)

      return
      end
      
      function wtsum(wt)
      implicit none
      integer p
      double precision wtt,wtsum
      complex*16 wt(10)

      include 'polarization.f'

      wtt=0d0

      do p=1,pol
         wtt=wtt+cdabs(wt(p))**2
      enddo

      wtsum=wtt

      return
      end

ccc   integrates bare + screened amplitude over k_t
ccc   (two-photon induced processes)
      subroutine schimcgamion_s2var(s2_in,p1x,p1y,p2x,p2y,outg)
      implicit none
      double precision x00p,wt,x00p2
      double precision xggmin,yp,ypmax,ypmin
      complex*16 zout,zout1,zoutg,zoutoff2,zoutoff2_nos2
      double precision tpx,tpy,tp2,t12,t22,t11,phiq
      double precision sc,qtmax,qt,screeningionint,qtmin
      double precision p1xp,p2xp,p1yp,p2yp
      double precision hy,hqt,hphi,del,dbl
      double precision p1x,p1y,p2x,p2y,denom,wtt_nos2,wtt_s2
      integer p,i,nphi,jqt,jphi
      complex*16 out(4,10),x0(10),outg(10),outg2(10)
      complex*16 zouttest,zouttest1,zout2,zoutoff2_2
      integer nk
      integer s2_in
      integer itvar
      common/itv/itvar

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
      include 'ion_inel.f'

      nphi=s2_in
      nk=s2_in

      wtt_s2=0d0
      wtt_nos2=0d0

      del=0.1d0

      qtmin=1d-7   
      qtmax=0.5d0

      if(ionbreakup)then
         if(int_01)then
            xggmin=1d-5         ! 01
            qtmax=0.1d0
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
         outg2(p)=0d0
         do i=1,4
            out(i,p)=0d0
         enddo
      enddo

      call wtgengam

      if(sfac)then

         do jqt=1,nk

            if(s2_in.eq.16)yp=((ypmax-ypmin)*x16(jqt)+ypmax+ypmin)/2d0
            if(s2_in.eq.32)yp=((ypmax-ypmin)*x32(jqt)+ypmax+ypmin)/2d0
            if(s2_in.eq.64)yp=((ypmax-ypmin)*x64(jqt)+ypmax+ypmin)/2d0
            if(s2_in.eq.96)yp=((ypmax-ypmin)*x96(jqt)+ypmax+ypmin)/2d0
            if(s2_in.eq.128)yp=((ypmax-ypmin)*x128(jqt)+ypmax+ypmin)/2d0
            if(s2_in.eq.256)yp=((ypmax-ypmin)*x256(jqt)+ypmax+ypmin)/2d0
            if(s2_in.eq.512)yp=((ypmax-ypmin)*x512(jqt)+ypmax+ypmin)/2d0
            if(s2_in.eq.1024)yp=((ypmax-ypmin)*x1024(jqt)+ypmax+ypmin)/2d0

            qt=dexp(yp)-xggmin**2
            
            tp2=qt
            qt=dsqrt(dabs(qt))

            sc=screeningionint(qt)


            zouttest=0d0
            zouttest1=0d0

            do jphi=1,nphi

               if(s2_in.eq.16)phiq=pi*(x16(jphi)+1d0)
               if(s2_in.eq.32)phiq=pi*(x32(jphi)+1d0)
               if(s2_in.eq.64)phiq=pi*(x64(jphi)+1d0)
               if(s2_in.eq.96)phiq=pi*(x96(jphi)+1d0)
               if(s2_in.eq.128)phiq=pi*(x128(jphi)+1d0)
               if(s2_in.eq.256)phiq=pi*(x256(jphi)+1d0)
               if(s2_in.eq.512)phiq=pi*(x512(jphi)+1d0)
               if(s2_in.eq.1024)phiq=pi*(x1024(jphi)+1d0)

               tpx=qt*dcos(phiq)
               tpy=qt*dsin(phiq)

               wt=(tp2+xggmin**2)/2d0*(ypmax-ypmin)/2d0*pi

               if(s2_in.eq.16)wt=wt*w16(jphi)*w16(jqt)
               if(s2_in.eq.32)wt=wt*w32(jphi)*w32(jqt)
               if(s2_in.eq.64)wt=wt*w64(jphi)*w64(jqt)
               if(s2_in.eq.96)wt=wt*w96(jphi)*w96(jqt)
               if(s2_in.eq.128)wt=wt*w128(jphi)*w128(jqt)
               if(s2_in.eq.256)wt=wt*w256(jphi)*w256(jqt)
               if(s2_in.eq.512)wt=wt*w512(jphi)*w512(jqt)
               if(s2_in.eq.1024)wt=wt*w1024(jphi)*w1024(jqt)

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
                  if(ion_inel)then
                  call formfacgamion_inel(1,t12,t22,x00p)
                  call formfacgamion_inel(2,t12,t22,x00p2)
                  else
                  call formfacgamionp(1,t12,t22,x00p)
                  call formfacgamionp(2,t12,t22,x00p2)
                  endif
               endif

           do p=1,pol  

              if(offshell)then

               if(ion_inel)then
         call formfacgamoff_ionp_surv(1,p,p1xp,p1yp,p2xp,p2yp,zout1)
         call formfacgamoff_ionp_surv(2,p,p1xp,p1yp,p2xp,p2yp,zout2)
         outg2(p)=outg2(p)+wt*sc*zout2

               else
               call formfacgamoff_ion_surv(p,p1xp,p1yp,p2xp,p2yp,zout1)
               endif

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
         if(ion_inel)then
         if(ion_em)then
         call formfacgamion_ion_em(1,t11,t22,x00p)
         call formfacgamion_ion_em(2,t11,t22,x00p2)
         else
         call formfacgamion_inel(1,t11,t22,x00p)
         call formfacgamion_inel(2,t11,t22,x00p2)
         endif
         else
         call formfacgamionp(1,t11,t22,x00p)
         call formfacgamionp(2,t11,t22,x00p2)
         endif
      endif

      do p=1,pol  

         if(ion_em)then


         call formfacgamoff(p,p1x,p1y,p2x,p2y,zoutoff2)
         outg(p)=dsqrt(cdabs(zoutoff2)**2)
      
         elseif(offshell)then

            if(ion_inel)then
            if(sfac)then
            call formfacgamoff_ionp_surv(1,p,p1x,p1y,p2x,p2y,zoutoff2)
            call formfacgamoff_ionp_surv(2,p,p1x,p1y,p2x,p2y,zoutoff2_2)
            call formfacgamoff_ionp(p,p1x,p1y,p2x,p2y,zoutoff2_nos2)
            else
            call formfacgamoff_ionp(p,p1x,p1y,p2x,p2y,zoutoff2)
            endif
            else
            if(sfac)then
            call formfacgamoff_ion_surv(p,p1x,p1y,p2x,p2y,zoutoff2)
            call formfacgamoff_ion(p,p1x,p1y,p2x,p2y,zoutoff2_nos2)
            else
            call formfacgamoff_ion(p,p1x,p1y,p2x,p2y,zoutoff2)
            endif
            endif

            if(ionbreakup)then
               if(sfac)then
               if(wrho)then
               else
                  if(fAA.eq.'00'.or.fAA.eq.'AA'.or.fAA.eq.'A0')then ! exp(-omega)*p and not 1-exp(-omega)*p (more stable) for other cases
                     outg(p)=outg(p)+zoutoff2
                     if(cdabs(zoutoff2).lt.cdabs(outg(p))*del)
     &                  outg(p)=zoutoff2
                     outg2(p)=outg2(p)+zoutoff2_2
                  endif
               endif
               else
                  outg(p)=outg(p)+zoutoff2
                  outg2(p)=outg2(p)+zoutoff2_2
               endif
            else

               outg(p)=outg(p)+zoutoff2 

               if(cdabs(zoutoff2).lt.cdabs(outg(p))*del)
     &           outg(p)=zoutoff2

               outg2(p)=outg2(p)+zoutoff2_2
            endif

            outg(p)=dsqrt(cdabs(outg(p))**2+cdabs(outg2(p)))
            if(sfac)then

            denom=dsqrt(cdabs(zoutoff2)**2+cdabs(zoutoff2_2)**2)

            if(denom.eq.0d0)then
            outg(p)=0d0

            else

            outg(p)=outg(p)*cdabs(zoutoff2_nos2)/denom

            endif
            endif

         else

            zout=-0.5d0*(ppa(p)+mma(p))*(p1x*p2x+p1y*p2y)
     &           -0.5d0*zi*(ppa(p)-mma(p))*(p1x*p2y-p2x*p1y)
     &           +0.5d0*(p1x*p2x-p1y*p2y+zi*(p1x*p2y+p1y*p2x))*mpa(p)
     &           +0.5d0*(p1x*p2x-p1y*p2y-zi*(p1x*p2y+p1y*p2x))*pma(p)

            if(ionbreakup)then
               if(sfac)then
               if(wrho)then
               else
                  if(fAA.eq.'00'.or.fAA.eq.'AA'.or.fAA.eq.'A0')then ! exp(-omega)*p and not 1-exp(-omega)*p (more stable) for other cases
                     outg(p)=outg(p)+zout*x00p*2d0
                  endif
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
