ccc   integrates bare + screened amplitude over k_t
ccc   (two-photon induced processes)
      subroutine schimcgam(p1x,p1y,p2x,p2y,outgg)
      implicit none
      complex*16 zoutoff2,zoutoff1,zoutg,zout1,zout2,zout
      double precision x00p2a,x00p2,x00p,sc,sc1
      double precision s2sd1,s2sd2,s2sd,s2out,s2i,s2dd
      double precision qtmax,qtcut,qt1,qt2,qt
      double precision phiq,hqt,hphi
      double precision p1xp,p2xp,p1yp,p2yp
      double precision del,dbl
      double precision p1x,p1y,p2x,p2y
      double precision wt,tpx,tpy,tp2,t11,t12,t22
      integer i1,i2,p,i,nphi,nqt,jqt,jphi
      complex*16 out(4,10),x0(10),outgg(10)
      complex*16 screen(2,2)
      complex*16 outg(2,10)
      complex*16 outgb(10),outgbo(10)
      complex*16 zouts_1arr(10),zouts_2arr(10)

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
      include 'beam.f'
      include 'mom.f'
      include 'diss.f'
      include 'gamma.f'
      include 'x.f'
      include 'diff.f'
      include 'gaussvars.f'
      include 'xb.f'
      include 'mp.f'

c      do p=1,pol
c         outgg(p)=1d0
c      enddo
c      return

      del=0.5d0

      nphi=s2int
      nqt=s2int*2

      qtmax=2d0
      qtcut=1d0

      hqt=qtmax/dble(nqt)
      hphi=2d0*pi/dble(nphi)

      zoutg=0d0

      do p=1,pol
         outgg(p)=0d0
         outg(1,p)=0d0
         outg(2,p)=0d0
         do i=1,4
            out(i,p)=0d0
         enddo
      enddo

      call wtgengam

      if(sfac)then

         if(offshell)then
            t11=p1x**2+p1y**2
            t22=p2x**2+p2y**2
            qt1=dsqrt(t11)
            qt2=dsqrt(t22)


            if(diff.eq.'sd')then
               if(diss1)then
                  if(qt1.gt.qtcut)goto 777
               elseif(diss2)then
                  if(qt2.gt.qtcut)goto 777
               endif
            elseif(diff.eq.'dd')then
               if(qt1.gt.qtcut)goto 777
               if(qt2.gt.qtcut)goto 777
            endif
         endif

         do jqt=1,nqt

            qt=qtmax/2d0*(xikt(jqt)+1d0)

            tp2=qt**2

            do i1=1,nch
               do i2=1,nch
                  call screeningint(i1,i2,tp2,sc,sc1)
                  screen(i1,i2)=sc
               enddo
            enddo

            do jphi=1,nphi

               phiq=pi*(xiphib(jphi)+1d0)

               tpx=qt*dcos(phiq)
               tpy=qt*dsin(phiq)

               wt=qt
               wt=wt*wiphib(jphi)*wikt(jqt)*pi*qtmax/2d0

               p1xp=p1x-tpx
               p1yp=p1y-tpy
               t12=p1xp**2+p1yp**2
               p2xp=tpx+p2x
               p2yp=tpy+p2y
               t22=p2xp**2+p2yp**2

               if(offshell.eqv..false.)then
                  call formfacgam(1,t12,t22,x00p)
                  call formfacgam(2,t12,t22,x00p2a)
               endif

               do p=1,pol

                  if(offshell)then

                     if(proc.eq.54.or.proc.eq.55)then
                        if(p.eq.4)then
                           zout1=zouts_1arr(1)
                           zout2=zouts_2arr(1)
                           goto 111
                        endif
                        if(p.eq.3)then
                           zout1=zouts_1arr(2)
                           zout2=zouts_2arr(2)
                           goto 111
                        endif
                        if(p.eq.6)then
                           zout1=zouts_1arr(5)
                           zout2=zouts_2arr(5)
                           goto 111
                        endif
                        if(p.eq.8)then
                           zout1=zouts_1arr(7)
                           zout2=zouts_2arr(7)
                           goto 111
                        endif
                     endif


                     call formfacgamoff_surv(1,p,p1xp,p1yp,p2xp,p2yp
     &                    ,zout1)
                     call formfacgamoff_surv(2,p,p1xp,p1yp,p2xp,p2yp
     &                    ,zout2)

                     zouts_1arr(p)=zout1
                     zouts_2arr(p)=zout2

 111                 x00p2a=dble(zout2)
                  endif

              do i1=1,nch
                 do i2=1,nch

                    x0(p)=pp0(i1)*pp0(i2)/dble(nch)**2/gaa(i1)
     &                   /gaa(i2)

                    if(offshell)then
                       zout=zout1
                    else
                       zout=-0.5d0*(ppa(p)+mma(p))*(p1xp*p2xp+p1yp*p2yp)
     &                   -0.5d0*zi*(ppa(p)-mma(p))*(p1xp*p2yp-p2xp*p1yp)
     &                      +0.5d0*(p1xp*p2xp-p1yp*p2yp
     &                      +zi*(p1xp*p2yp+p1yp*p2xp))*mpa(p)
     &                      +0.5d0*(p1xp*p2xp-p1yp*p2yp
     &                      -zi*(p1xp*p2yp+p1yp*p2xp))*pma(p)

                       zout=zout*2d0
                       zout=zout*x00p

                    endif

                       outg(1,p)=outg(1,p)+x0(p)*wt
     &                      *screen(i1,i2)*zout

                       if(offshell)then

                          outg(2,p)=outg(2,p)+x0(p)*wt
     &                      *screen(i1,i2)*zout2

                       elseif(p.eq.1)then
                          x00p2=x00p2a*pp0(i1)*pp0(i2)/dble(nch)**2
     &                         /gaa(i1)/gaa(i2)

                          zoutg=zoutg+x00p2*wt*screen(i1,i2)

                       endif

                 enddo
              enddo
           enddo
      enddo

      enddo


      endif


 777  t11=p1x**2+p1y**2
      t22=p2x**2+p2y**2


      if(offshell.eqv..false.)then
      if(beam.eq.'prot')then
         call formfacgam(1,t11,t22,x00p)
         call formfacgam(2,t11,t22,x00p2)
      elseif(beam.eq.'el')then
         call formfacgamel(1,t11,t22,x00p)
         call formfacgamel(2,t11,t22,x00p2)
      endif
      endif

      do p=1,pol

         if(offshell)then
c            call formfacgamoff(1,p,p1x,p1y,p2x,p2y,zoutoff1)
            call formfacgamoff(p,p1x,p1y,p2x,p2y,zoutoff2)

            outgb(p)=dsqrt(cdabs(zoutoff2)**2)
c            outgb(p)=dsqrt(cdabs(zoutoff1)**2+cdabs(zoutoff2)**2)

            if(sfac)then

               if(proc.eq.54.or.proc.eq.55)then
                  if(p.eq.4)then
                     zoutoff1=zouts_1arr(1)
                     zoutoff2=zouts_2arr(1)
                     goto 112
                  endif
                  if(p.eq.3)then
                     zoutoff1=zouts_1arr(2)
                     zoutoff2=zouts_2arr(2)
                     goto 112
                  endif
                  if(p.eq.6)then
                     zoutoff1=zouts_1arr(5)
                     zoutoff2=zouts_2arr(5)
                     goto 112
                  endif
                  if(p.eq.8)then
                     zoutoff1=zouts_1arr(7)
                     zoutoff2=zouts_2arr(7)
                     goto 112
                  endif
               endif

            call formfacgamoff_surv(1,p,p1x,p1y,p2x,p2y,zoutoff1)
            call formfacgamoff_surv(2,p,p1x,p1y,p2x,p2y,zoutoff2)


            zouts_1arr(p)=zoutoff1
            zouts_2arr(p)=zoutoff2

 112        outgbo(p)=dsqrt(cdabs(zoutoff1)**2+cdabs(zoutoff2)**2)

            outg(1,p)=outg(1,p)+zoutoff1
            outg(2,p)=outg(2,p)+zoutoff2

            if(cdabs(zoutoff1).lt.cdabs(outg(1,p))*del)
     &           outg(1,p)=zoutoff1
            if(cdabs(zoutoff2).lt.cdabs(outg(2,p))*del)
     &           outg(2,p)=zoutoff2

            dbl=dsqrt(cdabs(outg(1,p))**2+cdabs(outg(2,p))**2)
            outgg(p)=dbl

            endif

         else
            zout=-0.5d0*(ppa(p)+mma(p))*(p1x*p2x+p1y*p2y)
     &           -0.5d0*zi*(ppa(p)-mma(p))*(p1x*p2y-p2x*p1y)
     &           +0.5d0*(p1x*p2x-p1y*p2y+zi*(p1x*p2y+p1y*p2x))*mpa(p)
     &           +0.5d0*(p1x*p2x-p1y*p2y-zi*(p1x*p2y+p1y*p2x))*pma(p)

            outg(1,p)=outg(1,p)+zout*x00p*2d0
            dbl=cdabs(outg(1,p))

            dbl=dsqrt(dbl**2+cdabs(zoutg+x00p2)**2*pincarr(p))
            outgb(p)=dbl
            outgg(p)=dbl
         endif

      enddo

      if(sfac.eqv..false.)then
         do p=1,pol
            outgg(p)=outgb(p)
         enddo
         return
      endif

      if(offshell)then

         qt1=dsqrt(t11)
         qt2=dsqrt(t22)

         if(diff.eq.'el')then

            do p=1,pol
               if(cdabs(outgbo(p)).lt.1d-60)then
                  s2i=1d0
               else
                  s2i=cdabs(outgg(p)/outgbo(p))**2
               endif

               if(sfac)then
                  outgg(p)=outgb(p)*dsqrt(s2i)
               else
                  outgg(p)=outgb(p)
               endif

            enddo

         elseif(diff.eq.'sd')then

            do p=1,pol
               if(cdabs(outgbo(p)).lt.1d-60)then
                  s2i=1d0
               else
                  s2i=cdabs(outgg(p)/outgbo(p))**2
               endif

               if(diss1)then
                  call s2_diss_int(1,qt2,x2,s2sd)
                  call s2exp(qt1,s2i,s2sd,s2out)
               elseif(diss2)then
                  call s2_diss_int(1,qt1,x1,s2sd)
                  call s2exp(qt2,s2i,s2sd,s2out)
               endif

               if(sfac)then
                  outgg(p)=outgb(p)*dsqrt(s2out)
               else
                  outgg(p)=outgb(p)
               endif

            enddo

         elseif(diff.eq.'dd')then

            do p=1,pol

               if(cdabs(outgbo(p)).lt.1d-60)then
                  s2i=1d0
               else
                  s2i=cdabs(outgg(p)/outgbo(p))**2
               endif

               call s2_diss_int(2,qt1,x1,s2sd1)
               call s2_diss_int(2,qt2,x2,s2sd2)
               call s2_diss_int(3,qt2,x2,s2dd)

               call s2exp2(qt1,qt2,s2i,s2sd1,s2sd2,s2dd,s2out)

               if(sfac)then
                  outgg(p)=outgb(p)*dsqrt(s2out)
               else
                  outgg(p)=outgb(p)
               endif

            enddo

         endif


      endif


      return
      end
