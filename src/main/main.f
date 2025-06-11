c   calculates CEP cross section
      function cs(rarr,wgt)
      implicit none
      complex*16 wt(10),wtn(10),wtd(10),wtpvar(3,10),wt_0(10)
     &,wt_atau(10),wt_atauonly_lin(10)
      double precision rarr(10),wtr(10)
      integer i,p,icut
      double precision rphi,ran2,cs
      double precision ypp,ypmin1,ypmin2,ypmin,ypmax,ypmax1,ypmax2,
     &     yp,ymin1,ymax1,ycut,xgmin,wty1,wty2,wty,wtt,wtpt,wtpol,wtt_0,
     &     wtt_atau,wtt_atauonly_lin
      double precision wtdiss1,wtdiss2,wtc0,wt4,wt6,wt3a,wt3b,wt3,wt2a,
     &     wt2b,wt2,wt1
      double precision xglu,wpsi
      double precision val,ry
      double precision root1sq,root2sq,rmx,rm,cc1,cc2,aa1,aa2
      double precision rdisst,rdiss1,rdiss2,rdiss,r1,r2,r3,r4,r5
      double precision qsq1,qsq2,qsq1tt,qsq2tt
      double precision ptxx,ptxsq,ptmin,ptmax1,ptmax2,ptmax,ptdif
      double precision pt2x,pt2y,pt1y,pt1x,pt1sq,pt2sq,phi2,phi1
      double precision ps,p2p,p2m,p1p,p1m,ktcut
      double precision msub,mpp1,mpp2,mdissmax,lmdissmin,lmdissmax,
     &     lmdiss1,lmdiss2
      double precision jrho,jmono,jchi,jalp
      double precision wgt

      include 'polvecs.f'
      include 'gencuts.f'
      include 'pi.f'
      include 'unweighted.f'
      include 'bin.f'
      include 'x.f'
      include 'zi.f'
      include 'vars.f'
      include 'survpars.f'
      include 'mom.f'
      include 'mandelstam.f'
      include 'range.f'
      include 'polarization.f'
      include 'proc.f'
      include 'mq.f'
      include 'norm.f'
      include 'meta.f'
      include 'mres.f'
      include 'decay.f'
      include 'mp.f'
      include 'quarkonia.f'
      include 'scorr.f'
      include 'photo.f'
      include 'bpsi.f'
      include 'mpip.f'
      include 'mkp.f'
      include 'gamma.f'
      include 'ewpars.f'
      include 'eff.f'
      include 'brs.f'
      include 'wmax.f'
      include 'unwsurv.f'
      include 'survin.f'
      include 'widths.f'
      include 'surv.f'
      include 'beam.f'
      include 'mion.f'
      include 'ionqcd.f'
      include 's2qcd.f'
      include 'ion.f'
      include 'rho0.f'
      include 'mcharg.f'
      include 'xb.f'
      include 'diss.f'
      include 'diff.f'
      include 'pdg.f'
      include 'polww.f'
      include 'enew.f'
      include 'wwpars.f'
      include 'elcollw.f'
      include 'phimu.f'
      include 'wdecay.f'
      include 'p0Xn.f'
      include 'mxs.f'
      include 'tau.f'
      
      wtt=0d0

      elcollw=.false.

      if(beam.eq.'ionp')call pAinit

      if(dps.eq.1)then
         mx=mres
         if(proc.eq.48)call bwrho(mx,jrho)
         if(fwidth)then
            if(proc.gt.20.and.proc.lt.38)then
               call bwchi(mx,jchi)
            endif
            if(proc.eq.68)call bwalp(mx,jalp)
            if(proc.eq.69.or.proc.eq.70)then
               rm=rarr(5)
               mx=mmin**2+(mmax**2-mmin**2)*rm
               mx=dsqrt(mx)
               call bwmono(mx,jmono)
               jmono=jmono*(mmax**2-mmin**2)
            endif
         endif
      else
         if(photo.or.gamma)then
         else
            if(mmin.lt.2d0)then
               print*,'WARNING: need mmin > 2 GeV for QCD processes'
               print*,'Please edit in input card'
               STOP 1
            endif
         endif
         if(mmax.gt.rts)mmax=rts
         if(proc.eq.82)then
            if(mmin.lt.mxs)mmin=mxs
         elseif(proc.eq.83.or.proc.eq.84)then
            if(mmin.lt.mxs+mq)mmin=mxs+mq
         else
            if(mmin.lt.2d0*mq)mmin=2d0*mq
         endif
         rm=rarr(5)
         mx=mmin+(mmax-mmin)*rm
         if(mmin.lt.1d-3)mmin=1d-3
         msub=1d0/1d0*(1d0/mmax**1d0+
     &           (1d0/mmin**1d0-1d0/mmax**1d0)*rm)
         mx=(1d0/msub/1d0)**(1d0/1d0)
      endif

      if(beam.eq.'prot')then
         mpp1=mp
         mpp2=mp
      elseif(beam.eq.'ion')then
         mpp1=mion
         mpp2=mion
      elseif(beam.eq.'ionp')then
         mpp1=mp
         mpp2=mion
      elseif(beam.eq.'el')then
         mpp1=me
         mpp2=me
      endif

      if(beam.eq.'prot'.and.difftot)then
         if(dps.eq.2)then
            rdisst=rarr(9)
         else
            rdisst=rarr(8)
         endif
         if(rdisst.lt.1d0/3d0)then
            diff='el'
            diss1=.false.
            diss2=.false.
            elcoll=.false.
         elseif(rdisst.lt.2d0/3d0)then
            diff='sd'
            if(unw)elcoll=.true.
         else
            diff='dd'
            diss1=.true.
            diss2=.true.
            elcoll=.false.
         endif
      endif

      if(photo)then

         r1=rarr(2)
         r3=rarr(3)
         r4=rarr(4)
         r5=rarr(5)

         r2=ran2()

         ptmax=dsqrt(5d0)
         ptmin=0d0

         pt2sq=(ptmax-ptmin)*r1+ptmin
         pt2sq=pt2sq**2

         phi1=2d0*pi*r2
         phi2=2d0*pi*r3+phi1

         rmx=dsqrt(pt2sq+mx**2)
         wtpt=2d0*dsqrt(pt2sq)*(ptmax-ptmin)
         xgmin=(mx/rts)**2

         if(r5.gt.0.5d0)then
            ypmax=dlog(xgmin**2*mpp1**2+ptmax**2)
            ypmin=dlog(xgmin**2*mpp1**2)
            yp=(ypmax-ypmin)*r4+ypmin
            pt1sq=dexp(yp)-xgmin**2*mpp1**2
            wt1=xgmin**2*mpp1**2+pt1sq
         else
            ypmax=dlog(xgmin**2*mpp2**2+ptmax**2)
            ypmin=dlog(xgmin**2*mpp2**2)
            yp=(ypmax-ypmin)*r4+ypmin
            pt1sq=pt2sq
            pt2sq=dexp(yp)-xgmin**2*mpp2**2
            wt1=xgmin**2*mpp2**2+pt2sq
         endif

         pt2x=dsqrt(pt2sq)*dcos(phi2)
         pt2y=dsqrt(pt2sq)*dsin(phi2)
         pt1x=dsqrt(pt1sq)*dcos(phi1)
         pt1y=dsqrt(pt1sq)*dsin(phi1)

         ptxx=(pt1x+pt2x)**2+(pt1y+pt2y)**2
         rmx=dsqrt(ptxx+mx**2)

         if(beam.eq.'ionp')then
            ymax=dlog(rts/rmx)
            ymin=-ymax
         endif

         ry=rarr(1)
         yx=ymin+(ymax-ymin)*ry

         wty=ymax-ymin

         if(beam.eq.'ionp')r5=0d0

         if(r5.gt.0.5d0)then    ! photon emitted from q(1,k)
            xgam=rmx*dexp(yx)/rts ! photon mom. fraction
            wpsi=dsqrt(xgam*s)  ! proton-photon cms energy
            xglu=(rmx/wpsi)**2  ! gluon mom. fraction
            x1=xgam
            x2=xglu
            prot=1
         else                   ! photon emitted from q(1,k)
            xgam=rmx*dexp(-yx)/rts
            wpsi=dsqrt(xgam*s)
            xglu=(rmx/wpsi)**2
            x2=xgam
            x1=xglu
            prot=2
         endif

         bpsi=bpsi0+4d0*alphapb*dlog(wpsi/w0b)

         if(x1.gt.1d0)goto 777
         if(x2.gt.1d0)goto 777

      elseif(gamma)then


         if(beam.eq.'prot'.and.diff.eq.'sd')then
            rdiss=ran2()
            if(rdiss.gt.0.5d0)then
               diss1=.true.
               diss2=.false.
            else
               diss1=.false.
               diss2=.true.
            endif
            if(diffsd.eq.'sda')then
               diss1=.true.
               diss2=.false.
            elseif(diffsd.eq.'sdb')then
               diss1=.false.
               diss2=.true.
            endif

         endif

         r2=rarr(2)
         r3=rarr(3)
         r4=rarr(4)

         r1=ran2()


         phi1=2d0*pi*r1
         phi2=2d0*pi*r2+phi1

         if(beam.eq.'prot'.or.beam.eq.'ion'.or.beam.eq.'ionp')then
            ptmax=dsqrt(10d0)
         elseif(beam.eq.'el')then
            ptmax=dsqrt(50d0)
         endif

         ptmin=0d0
         xgmin=(mx/rts)**2

         if(diss1)then
            ptmax1=rts/2d0
            ypmax1=dlog(xgmin**2*mpp1**2+ptmax1**2)
            ypmin1=dlog(xgmin**2*mpp1**2)
            yp=(ypmax1-ypmin1)*r3+ypmin1
            pt1sq=dexp(yp)-xgmin**2*mpp1**2
            wty1=(ypmax1-ypmin1)*(xgmin**2*mpp1**2+pt1sq)
         else
            ypmax1=dlog(xgmin**2*mpp1**2+ptmax**2)
            ypmin1=dlog(xgmin**2*mpp1**2)
            yp=(ypmax1-ypmin1)*r3+ypmin1
            pt1sq=dexp(yp)-xgmin**2*mpp1**2
            wty1=(ypmax1-ypmin1)*(xgmin**2*mpp1**2+pt1sq)
         endif

         if(diss2)then
            ptmax2=rts/2d0
            ypmax2=dlog(xgmin**2*mpp2**2+ptmax2**2)
            ypmin2=dlog(xgmin**2*mpp2**2)
            ypp=(ypmax2-ypmin2)*r4+ypmin2
            pt2sq=dexp(ypp)-xgmin**2*mpp2**2
            wty2=(ypmax2-ypmin2)*(xgmin**2*mpp2**2+pt2sq)
         else
            ypmax2=dlog(xgmin**2*mpp2**2+ptmax**2)
            ypmin2=dlog(xgmin**2*mpp2**2)
            ypp=(ypmax2-ypmin2)*r4+ypmin2
            pt2sq=dexp(ypp)-xgmin**2*mpp2**2
            wty2=(ypmax2-ypmin2)*(xgmin**2*mpp2**2+pt2sq)
         endif

ccccc

         if(pt1sq.lt.0d0)then
            pt1sq=0d0
         endif

         if(pt2sq.lt.0d0)then
            pt2sq=0d0
         endif

         pt2x=dsqrt(pt2sq)*dcos(phi2)
         pt2y=dsqrt(pt2sq)*dsin(phi2)
         pt1x=dsqrt(pt1sq)*dcos(phi1)
         pt1y=dsqrt(pt1sq)*dsin(phi1)

         ptxx=(pt1x+pt2x)**2+(pt1y+pt2y)**2
         rmx=dsqrt(ptxx+mx**2)


         if(beam.eq.'ionp')then
            ymax=dlog(rts/rmx)
            ymin=-ymax
         endif

         ry=rarr(1)
         yx=ymin+(ymax-ymin)*ry

         wty=ymax-ymin

         x1=rmx*dexp(yx)/rts    ! photon 1 mom. fraction
         x2=rmx*dexp(-yx)/rts   ! photon 2 mom. fraction

         if(x1.gt.1d0)goto 777
         if(x2.gt.1d0)goto 777

         wpsi=dsqrt(x1*x2*s)

      else

         r1=rarr(2)
         r2=rarr(3)
         r4=rarr(4)

         r3=ran2()

         ptmax=dsqrt(2d0)
         ptmin=0d0

         pt1sq=r1*ptmax
         pt2sq=r2*ptmax

         pt1sq=ptmin+(ptmax-ptmin)*r1
         pt2sq=ptmin+(ptmax-ptmin)*r2

         pt1sq=pt1sq**2
         pt2sq=pt2sq**2

         phi1=2d0*pi*r3
         phi2=2d0*pi*r4+phi1

         pt1x=dsqrt(pt1sq)*dcos(phi1)
         pt1y=dsqrt(pt1sq)*dsin(phi1)
         pt2x=dsqrt(pt2sq)*dcos(phi2)
         pt2y=dsqrt(pt2sq)*dsin(phi2)

         ptxsq=(pt1x+pt2x)**2+(pt1y+pt2y)**2
         rmx=dsqrt(mx**2+ptxsq)

cccccccccccccccccccccc

         ymax1=ymax
         ymin1=ymin
         ycut=dlog(rts/rmx)
         ry=rarr(1)
         yx=ymin+(ymax-ymin)*ry

         ymin=ymin1
         ymax=ymax1

         x1=rmx/rts*dexp(yx)
         x2=rmx/rts*dexp(-yx)

         if(x1.gt.1d0)goto 777
         if(x2.gt.1d0)goto 777

      endif

cccccccc   Dissociation



      if(beam.eq.'prot')then

         mdissmax=rts

         lmdissmax=dlog(mdissmax)
         lmdissmin=dlog(mp)

      if(diss1)then
         if(dps.eq.2)then
            rdiss1=rarr(7)
         else
            rdiss1=rarr(5)
         endif

         lmdiss1=lmdissmin+(lmdissmax-lmdissmin)*rdiss1
         mdiss1=dexp(lmdiss1)
         wtdiss1=2d0*mdiss1**2*(lmdissmax-lmdissmin)

      else
         mdiss1=mp
         wtdiss1=1d0
      endif

      if(diss2)then
         if(dps.eq.2)then
            rdiss2=rarr(8)
            if(diff.eq.'sd')rdiss2=rarr(7)
         else
            if(diff.eq.'sd')then
               rdiss2=rarr(5)
            else
               rdiss2=rarr(6)
            endif
         endif

         lmdiss2=lmdissmin+(lmdissmax-lmdissmin)*rdiss2
         mdiss2=dexp(lmdiss2)
         wtdiss2=2d0*mdiss2**2*(lmdissmax-lmdissmin)

      else
         mdiss2=mp
         wtdiss2=1d0
      endif

      mpp1=mdiss1
      mpp2=mdiss2

      else

         wtdiss1=1d0
         wtdiss2=1d0

      endif

cccccccccccc

      aa1=(1d0-x1)*rts/dsqrt(2d0)
      aa2=(1d0-x2)*rts/dsqrt(2d0)
      cc1=0.5d0*(pt2sq+mpp2**2)
      cc2=0.5d0*(pt1sq+mpp1**2)

c     impose massive on-shell condition by solving
c                   p1+ + cc1/p2- = aa1
c                   p2- + cc2/p1+ = aa2

      root1sq=(cc1-cc2-aa1*aa2)**2-4d0*cc2*aa1*aa2
      root2sq=(cc2-cc1-aa1*aa2)**2-4d0*cc1*aa1*aa2

      if(root1sq.le.0d0.or.root2sq.le.0d0)goto 777

      p1p=(cc2-cc1+aa1*aa2+dsqrt(root1sq))/(2d0*aa2)
      p2m=(cc1-cc2+aa1*aa2+dsqrt(root2sq))/(2d0*aa1)
      p1m=(pt1sq+mpp1**2)/(2d0*p1p)
      p2p=(pt2sq+mpp2**2)/(2d0*p2m)

      if(p1m.lt.0d0)goto 777
      if(p2p.lt.0d0)goto 777

      q(1,3)=pt1x
      q(2,3)=pt1y
      q(3,3)=(p1p-p1m)/dsqrt(2d0)
      q(4,3)=(p1p+p1m)/dsqrt(2d0)

      q(1,4)=pt2x
      q(2,4)=pt2y
      q(3,4)=(p2p-p2m)/dsqrt(2d0)
      q(4,4)=(p2p+p2m)/dsqrt(2d0)


      do i=1,4
         q(i,5)=q(i,1)+q(i,2)-q(i,3)-q(i,4)
      enddo

      if(beam.eq.'prot')then

         qsq1tt=(q(4,3)-q(4,1))**2-(q(3,3)-q(3,1))**2-(q(2,3)-q(2,1))**2
     &        -(q(1,3)-q(1,1))**2

         qsq2tt=(q(4,4)-q(4,2))**2-(q(3,4)-q(3,2))**2-(q(2,4)-q(2,2))**2
     &        -(q(1,4)-q(1,2))**2

         qsq1=-qsq1tt
         qsq2=-qsq2tt

ccccccc

         if(diss1)then
            xb1=qsq1/(qsq1+mdiss1**2-mp**2)
         else
            xb1=1d0
         endif
         if(diss2)then
            xb2=qsq2/(qsq2+mdiss2**2-mp**2)
         else
            xb2=1d0
         endif

         if(offshell)then
         else
            if(qsq1.gt.mx**2)goto 777
            if(qsq2.gt.mx**2)goto 777
         endif

         if(diss1.and.xb1.gt.1d0)goto 777
         if(diss2.and.xb2.gt.1d0)goto 777

         if(diss1)wtdiss1=wtdiss1/(qsq1+mdiss1**2-mp**2)
         if(diss2)wtdiss2=wtdiss2/(qsq2+mdiss2**2-mp**2)

      endif

ccccccccccccccccccccccccccccccccc

c      rphi=rarr(8)+r1+r2

      if(dps.eq.2)then
         rphi=ran2()
         if(proc.eq.82)then
            call twojetps_rm(mx,0d0,mxs,rarr(6),rphi,ps,uh,th)
         elseif(proc.eq.83.or.proc.eq.84)then
            call twojetps_rm(mx,mq,mxs,rarr(6),rphi,ps,uh,th)
         else
            call twojetps(mx,mq,rarr(6),rphi,ps,uh,th)
         endif
      elseif(dps.eq.3)then
         call threejetps(mx,mq,ps)
      elseif(dps.eq.12)then
         if(proc.eq.19)then
            call twojetpsm(mx,mpsi,mpsip,ps,uh,th)
         else
            call twojetpsm(mx,meta,metap,ps,uh,th)
         endif
      endif

cccccccccccccccccc

cccc decays

ccccccccccccccccccc

         wt2=1d0
         wt3=1d0
         wt4=1d0
         wt6=1d0

         if(proc.eq.1.or.proc.eq.60.or.proc.eq.67)then
            call twobody(1,5,6,7,mb,mb,wt2)
         elseif(proc.eq.68.or.proc.eq.69.or.proc.eq.70)then
            call twobody(1,5,6,7,0d0,0d0,wt2)
         elseif(proc.eq.18.or.proc.eq.19.or.proc.eq.20)then
            call twobody(1,6,8,9,mmu,mmu,wt2a)
            call twobody(1,7,10,11,mmu,mmu,wt2b)
            wt2=wt2a*wt2b
         elseif(proc.eq.21.or.proc.eq.22.or.proc.eq.23)then
            call twobody(1,5,6,7,0d0,mpsi,wt2a)
            call twobody(2,7,8,9,mmu,mmu,wt2b)
            wt2=wt2a*wt2b
         elseif(proc.gt.23.and.proc.lt.29)then
            call twobody(1,5,6,7,m2b,m2b,wt2)
         elseif(proc.eq.29.or.proc.eq.30.or.proc.eq.31)then
            call fourbody(mpip,mpip,wt4)
         elseif(proc.eq.32.or.proc.eq.33.or.proc.eq.34)then
            call fourbody(mpip,mkp,wt4)
         elseif(proc.eq.35.or.proc.eq.36.or.proc.eq.37)then
            call sixbody(mpip,wt6)
         elseif(proc.eq.39.or.proc.eq.40.or.proc.eq.41)then
            call twobody(1,5,6,7,0d0,mups,wt2a)
            call twobody(2,7,8,9,mmu,mmu,wt2b)
            wt2=wt2a*wt2b
         elseif(proc.gt.41.and.proc.lt.47)then
            call twobody(1,5,6,7,m2b,m2b,wt2)
         elseif(proc.eq.48)then
            call twobody(1,5,6,7,mpip,mpip,wt2)
         elseif(proc.eq.49)then
            if(phimu)then
               call twobody(1,5,6,7,mmu,mmu,wt2)
            else
               call twobody(1,5,6,7,mkp,mkp,wt2)
            endif
         elseif(proc.eq.50.or.proc.eq.51.or.proc.eq.52)then
            if(psie)then
               call twobody(1,5,6,7,me,me,wt2)
            else
               call twobody(1,5,6,7,mmu,mmu,wt2)
            endif
         elseif(proc.eq.53)then
            call threebody(1,5,6,7,8,mpsi,mpip,mpip,wt3)
            call twobody(1,6,9,10,mmu,mmu,wt2a)
            wt3=wt3*wt2a
         elseif(proc.eq.54.or.proc.eq.61)then
            call twobodyw(6,8,9,0d0,mmu)
            call twobodyw(7,10,11,0d0,mmu)
         elseif(proc.eq.55.or.proc.eq.62)then
            if(wlp_lep)then
               call wwmix
            elseif(wlm_lep)then
               call wwmix
            endif 
            if(wlp.eq.'mu')then
               call twobodyw(6,8,9,0d0,mmu)
            elseif(wlp.eq.'el')then
               call twobodyw(6,8,9,0d0,me)
            elseif(wlp.eq.'had')then
               call twobodyw(6,8,9,mu_quark,md_quark)
            else
               call twobodyw(6,8,9,0d0,mtau)
            endif
            if(wlm.eq.'mu')then
               call twobodyw(7,10,11,0d0,mmu)
            elseif(wlm.eq.'el')then
               call twobodyw(7,10,11,0d0,me)
            elseif(wlm.eq.'had')then
               call twobodyw(7,10,11,mu_quark,md_quark)
            else
               call twobodyw(7,10,11,0d0,mtau)
            endif
c            print*,dsqrt(q(4,11)**2-q(3,11)**2-q(2,11)**2-q(1,11)**2)
         elseif(proc.eq.73)then
            call threebody(1,6,8,9,10,mneut,mmu,0d0,wt3a)
            call threebody(1,7,11,12,13,mneut,mmu,0d0,wt3b)
            wt3=wt3a*wt3b
         elseif(proc.eq.74)then
            call threebody(1,6,8,9,10,mneut,0d0,0d0,wt3a)
            call threebody(1,7,11,12,13,mneut,0d0,0d0,wt3b)
            wt3=wt3a*wt3b
         elseif(proc.eq.75)then
            call threebody(1,6,8,9,10,mneut,mmu,0d0,wt3a)
            call threebody(2,7,11,12,13,mneut,0d0,0d0,wt3b)
            wt3=wt3a*wt3b
         elseif(proc.eq.76)then
            call twobody(1,6,8,9,mneut,mmu,wt2a)
            call twobody(1,7,10,11,mneut,mmu,wt2b)
            wt2=wt2a*wt2b
         elseif(proc.eq.83)then
            call twobody(1,6,8,9,mmu,mmu,wt2a)
            wt2=wt2a
         elseif(proc.eq.84)then
            call twobody(1,6,8,9,me,me,wt2a)
            wt2=wt2a
         endif

         if(beam.eq.'ionp')call pAboost

ccccccccccccccccccc  cuts ccccccccccccccccc

         neff0=neff0+1

         if(gencuts)then
            call cut(icut)
            if(icut.eq.0)goto 777
         endif

         neff=neff+1

ccccccccccccccccccccccccccccccccccccccccccc

          if(proc.eq.22.or.proc.eq.25.or.proc.eq.27.or.proc.eq.30
     &        .or.proc.eq.33.or.proc.eq.36)then
             call genpol1(5,echi1)
          elseif(proc.eq.23.or.proc.eq.26.or.proc.eq.28.or.proc.eq.
     &            31.or.proc.eq.34.or.proc.eq.37)then
             call genpol2
          elseif(proc.eq.54.or.proc.eq.55)then
             call genpol1(6,ewp)
             call genpol1(7,ewm)
                if(qsq1.gt.1d0.and.mdiss1.gt.dsqrt(3.5d0))then
                   wgauge='unitary'
                endif
                if(qsq2.gt.1d0.and.mdiss2.gt.dsqrt(3.5d0))then
                   wgauge='unitary'
                endif
             if(wgauge.eq.'axial')then
                call genpolaxial1(6,ewp)
                call genpolaxial1(7,ewm)
             endif
             wgauge='axial'
          endif

          if(proc.eq.40.or.proc.eq.43.or.proc.eq.45)then
             call genpol1(5,echi1)
          elseif(proc.eq.41.or.proc.eq.44.or.proc.eq.46)then
             call genpol2
          endif

ccccccccc


         if(photo)then
             if(beam.eq.'prot')then
                call schimcphot(pt1x,pt1y,pt2x,pt2y,wt)
             elseif(beam.eq.'ionp')then
                call schimcphotionp(pt1x,pt1y,pt2x,pt2y,wt)
             elseif(beam.eq.'ion')then
                print*,'Photoproduction not currently available for AA'
                STOP 1
             endif
         elseif(gamma)then
            if(beam.eq.'prot'.or.beam.eq.'el')then
               if(deltau)then
                  atau_temp=atau
                  dtau_temp=dtau
                  atau=0d0
                  dtau=0d0
                  GC_6506=-dtau/2d0
                  GC_6515=dsqrt(pi/1.325070D+02)/mtau*atau/2d0*zi
                  call schimcgam(pt1x,pt1y,pt2x,pt2y,wt_0)
                  atau=atau_temp
                  dtau=dtau_temp
                  GC_6506=-dtau/2d0
                  GC_6515=dsqrt(pi/1.325070D+02)/mtau*atau/2d0*zi
                  atau_only=.true.
                  call schimcgam(pt1x,pt1y,pt2x,pt2y,wt_atau)
                  atau_only=.false.
               endif
               if(deltau.and.atau_quad)then
                  atau_only=.true.
                  atau_lin=.true.
                  atau_quad=.false.
                  call schimcgam(pt1x,pt1y,pt2x,pt2y,wt_atauonly_lin)
                  atau_only=.false.
                  atau_lin=.false.
                  atau_quad=.true.
               endif    
               if(int_atauonly)then
                  atau_only=.true.
                  atau_lin=.true.
                  atau_quad=.false.
                  call schimcgam(pt1x,pt1y,pt2x,pt2y,wt_0)
                  atau_lin=.false.
                  atau_quad=.true.
                  call schimcgam(pt1x,pt1y,pt2x,pt2y,wt_atau)
                  atau_lin=.false.
                  atau_quad=.false.
               endif
               call schimcgam(pt1x,pt1y,pt2x,pt2y,wt)
            elseif(beam.eq.'ion'.or.beam.eq.'ionp')then
               if(pAAvar)then
                  do p=1,pol
                     do i=1,3
                        wtpvar(i,p)=0d0
                     enddo
                  enddo
               if(sfac)then
                  do ifaa=1,3
                     if(ifaa.eq.1)then ! inclusive                                                       
                        ionbreakup=.false.
                     endif
                     if(ifaa.eq.2)then
                        ionbreakup=.true.
                        faa='00'
                     endif
                     if(ifaa.eq.3)then
                        ionbreakup=.true.
                        faa='XX'
                     endif
                     call schimcgamion(pt1x,pt1y,pt2x,pt2y,wt)
                     do p=1,pol
                        wtpvar(ifaa,p)=wt(p)
                     enddo
                  enddo
                  do p=1,pol
                     wtr(p)=cdabs(wtpvar(1,p))**2
     &               -cdabs(wtpvar(2,p))**2
     &               -cdabs(wtpvar(3,p))**2
                     wtr(p)=wtr(p)/2d0
                  enddo
               else
                  call schimcgamion(pt1x,pt1y,pt2x,pt2y,wt)
               endif
            else
               if(deltau)then
                  atau_temp=atau
                  dtau_temp=dtau
                  atau=0d0
                  dtau=0d0
                  GC_6515=dsqrt(pi/1.325070D+02)/mtau*atau/2d0*zi
                  GC_6506=-dtau/2d0
                  call schimcgamion(pt1x,pt1y,pt2x,pt2y,wt_0)
                  atau=atau_temp
                  dtau=dtau_temp
                  GC_6515=dsqrt(pi/1.325070D+02)/mtau*atau/2d0*zi
                  GC_6506=-dtau/2d0
                  atau_only=.true.
                  call schimcgamion(pt1x,pt1y,pt2x,pt2y,wt_atau)
                  atau_only=.false.
               endif
               if(deltau.and.atau_quad)then
                  atau_only=.true.
                  atau_lin=.true.
                  atau_quad=.false.
                  call schimcgamion(pt1x,pt1y,pt2x,pt2y,wt_atauonly_lin)
                  atau_only=.false.
                  atau_lin=.false.
                  atau_quad=.true.
               endif    
               if(int_atauonly)then
                  atau_only=.true.
                  atau_lin=.true.
                  atau_quad=.false.
                  call schimcgamion(pt1x,pt1y,pt2x,pt2y,wt_0)
                  atau_lin=.false.
                  atau_quad=.true.
                  call schimcgamion(pt1x,pt1y,pt2x,pt2y,wt_atau)
                  atau_lin=.false.
                  atau_quad=.false.
               endif
               call schimcgamion(pt1x,pt1y,pt2x,pt2y,wt)
            endif
         endif
         else
             call wtgen
             if(beam.eq.'prot')then
                call schimc(pt1x,pt1y,pt2x,pt2y,wt)
             elseif(beam.eq.'ion'.or.beam.eq.'ionp')then
                if(ionqcd.eq.'incoh')then
                   if(sfac)then
                      call schimc(pt1x,pt1y,pt2x,pt2y,wt)
                      do p=1,pol
                         if(beam.eq.'ion')wt(p)=wt(p)*dsqrt(s2qcd)*an
                         if(beam.eq.'ionp')wt(p)=wt(p)*dsqrt(s2qcd*an)
                      enddo
                   else
                      call schimcion(pt1x,pt1y,pt2x,pt2y,wt)
                      do p=1,pol
                         if(beam.eq.'ion')wt(p)=wt(p)*dsqrt(s2qcd)*an
                         if(beam.eq.'ionp')wt(p)=wt(p)*dsqrt(s2qcd*an)
                      enddo
                   endif
                elseif(ionqcd.eq.'coh')then
                   ptdif=dsqrt((pt1x-pt2x)**2+(pt1y-pt2y)**2)
                   ktcut=ptdif
                   call schimcion(pt1x,pt1y,pt2x,pt2y,wt)
                   if(sfac)then
                      sfac=.false.
                      call schimc(pt1x,pt1y,pt2x,pt2y,wtd)
                      sfac=.true.
                      call schimc(pt1x,pt1y,pt2x,pt2y,wtn)
                      do p=1,pol
                         wt(p)=wt(p)*cdabs(wtn(p))/cdabs(wtd(p))
                      enddo
                   endif
                endif
             endif
          endif

          wtt=0d0

          if(paavar)then
             if(sfac)then
                do p=1,pol
                   wtt=wtt+wtr(p)
                enddo
             else
                do p=1,pol
                   wtt=wtt+cdabs(wt(p))**2
                enddo
             endif
          else
            wtt_0=0d0
            wtt_atau=0d0
            wtt_atauonly_lin=0d0
             do p=1,pol
                wtt=wtt+cdabs(wt(p))**2
                if(deltau)wtt_0=wtt_0+cdabs(wt_0(p))**2
                if(deltau)wtt_atau=wtt_atau+cdabs(wt_atau(p))**2
                if(int_atauonly)wtt_0=wtt_0+cdabs(wt_0(p))**2
                if(int_atauonly)wtt_atau=wtt_atau+cdabs(wt_atau(p))**2
                if(deltau.and.atau_quad)
     &wtt_atauonly_lin=wtt_atauonly_lin+cdabs(wt_atauonly_lin(p))**2
             enddo
             if(deltau)wtt=wtt-wtt_0-wtt_atau
             if(int_atauonly)wtt=wtt-wtt_0-wtt_atau
             if(deltau.and.atau_quad)wtt
     &=wtt+wtt_atauonly_lin
         endif
 

         wtpol=1d0

         if(scorr)then
            if(proc.eq.18.or.proc.eq.19.or.proc.eq.20)then
               call jpsidecay(wt,wtt)
            endif
            if(proc.eq.21.or.proc.eq.39)then
               call chic0decay3(wtc0)
               wtt=wtt*wtc0
            endif
            if(proc.eq.22.or.proc.eq.40)then
               do i=4,6
                  wt(i)=conjg(wt(i-3))
               enddo
               call chic1decay3(wt,wtt)
            endif
           if(proc.eq.23.or.proc.eq.41)then
              do i=6,10
                 wt(i)=conjg(wt(i-5))
              enddo
              call chic2decay3(wt,wtt)
           endif
           if(proc.eq.25.or.proc.eq.43)then
              do i=4,6
                 wt(i)=conjg(wt(i-3))
              enddo
              call chic1decay2s(wt,wtt)
           endif
           if(proc.eq.26.or.proc.eq.44)then
              do i=6,10
                 wt(i)=conjg(wt(i-5))
              enddo
              call chic2decay2s(wt,wtt)
           endif
           if(proc.eq.27.or.proc.eq.45)then
              do i=4,6
                 wt(i)=conjg(wt(i-3))
              enddo
              call chic1decay2f(wt,wtt)
           endif
           if(proc.eq.28.or.proc.eq.46)then
              do i=6,10
                 wt(i)=conjg(wt(i-5))
              enddo
              call chic2decay2f(wt,wtt)
           endif
           if(proc.eq.50.or.proc.eq.51.or.proc.eq.52)then
              call jpsidecayphot(wtpol)
           endif
           if(proc.eq.54.or.proc.eq.55)then
              call wwcorr(wt,wtt)
           endif
         endif

         wtt=wtt*wt2*wt3*wt4*wt6


         if(decays)then
            do i=1,nbr
               wtt=wtt*br(i)
            enddo
         endif



         if(photo)then
            wtt=wtt*wty
            wtt=wtt*wtpt
            wtt=wtt*(ypmax-ypmin)*wt1
            wtt=wtt*(wpsi/w0)**delta*normp*bpsi
            if(scorr)wtt=wtt*wtpol
            if(proc.eq.48)wtt=wtt*jrho
         elseif(gamma)then
            wtt=wtt*wty1*wty2
            wtt=wtt*s/2d0/dabs(p1p*p2m-p1m*p2p)
            if(offshell)then
               if(dps.eq.2)then
                  wtt=wtt/(128d0*pi**5) ! Final-state PS measure dGamma
                  wtt=wtt*pi**2 ! photon phi integral
                  wtt=wtt/2d0   ! def of sigma
                  wtt=wtt*mx**2 ! dm^2/m^2 -> dm^2
                  wtt=wtt/s**2  ! One from def of sigma, one from dx_i -> dM
               elseif(dps.eq.1)then
                  wtt=wtt/(2d0*pi)**3 ! Final-state PS measure dGamma
                  wtt=wtt*pi**2 ! photon phi integral
                  wtt=wtt/2d0   ! def of sigma
                  wtt=wtt/s**2  ! One from def of sigma, one from dx_i -> dM
                  wtt=wtt/pi*mx**4 ! To match with below
               endif
            endif
            wtt=wtt*wty
            wtt=wtt*2d0/mx
            if(dps.eq.1)wtt=wtt*pi/2d0/mx**3
            if(dps.eq.2)wtt=wtt*mx**2*(1d0/mmin**1d0-1d0/mmax**1d0)/1d0
            if(scorr)wtt=wtt*wtpol
            if(fwidth)then
               if(proc.eq.69.or.proc.eq.70)wtt=wtt*jmono
            endif
            wtt=wtt*wtdiss1*wtdiss2
            if(diff.eq.'sd')wtt=wtt*2d0
            if(diffsd.eq.'sda'.or.diffsd.eq.'sdb')wtt=wtt/2d0
            if(difftot)wtt=wtt*3d0
            if(fwidth)then
               if(proc.eq.68)then
                  wtt=wtt*jalp
               endif
            endif

         else
            wtt=wtt*(ymax-ymin)
            wtt=wtt*4d0*ptmax**2*dsqrt(pt1sq*pt2sq)*pi**2
            if(fwidth)then
               if(proc.gt.20.and.proc.lt.38)then
                  wtt=wtt*jchi
               elseif(proc.eq.68)then
                  wtt=wtt*jalp
               endif
            endif
         endif

         wtt=wtt/sym

         if(photo)goto 888
         if(gamma)goto 888

ccccccccccccc 1 body phase space

         if(dps.eq.1)then
            wtt=wtt/(16d0**2*pi**5)
         endif

         if(dps.gt.1)then
            wtt=wtt/(64d0*pi**2)*ps
            wtt=wtt*2d0*mx
            wtt=wtt/(16d0**2*pi**6)
            wtt=wtt*mx**2*(1d0/mmin**1d0-1d0/mmax**1d0)/1d0
         endif

         wtt=wtt*conv*surv

 888     if((elcoll.eqv..true.).and.(unw.eqv..true.))then

            if(diss1)then
            else
               pt1sq=0d0
               pt1x=0d0
               pt1y=0d0
            endif

            if(diss2)then
            else
               pt2sq=0d0
               pt2x=0d0
               pt2y=0d0
            endif

            ptxsq=(pt1x+pt2x)**2+(pt1y+pt2y)**2
            rmx=dsqrt(mx**2+ptxsq)

            x1=rmx/rts*dexp(yx)
            x2=rmx/rts*dexp(-yx)

            if(x1.gt.1d0)goto 777
            if(x2.gt.1d0)goto 777

ccccccccccccccc

            aa1=(1d0-x1)*rts/dsqrt(2d0)
            aa2=(1d0-x2)*rts/dsqrt(2d0)
            cc1=0.5d0*(pt2sq+mpp2**2)
            cc2=0.5d0*(pt1sq+mpp1**2)

c     impose massive on-shell condition by solving
c     p1+ + cc1/p2- = aa1
c                   p2- + cc2/p1+ = aa2

            root1sq=(cc1-cc2-aa1*aa2)**2-4d0*cc2*aa1*aa2
            root2sq=(cc2-cc1-aa1*aa2)**2-4d0*cc1*aa1*aa2

            if(root1sq.le.0d0.or.root2sq.le.0d0)goto 777

            p1p=(cc2-cc1+aa1*aa2+dsqrt(root1sq))/(2d0*aa2)
            p2m=(cc1-cc2+aa1*aa2+dsqrt(root2sq))/(2d0*aa1)
            p1m=(pt1sq+mpp1**2)/(2d0*p1p)
            p2p=(pt2sq+mpp2**2)/(2d0*p2m)

            if(p1m.lt.0d0)goto 777
            if(p2p.lt.0d0)goto 777

            q(1,3)=pt1x
            q(2,3)=pt1y
            q(3,3)=(p1p-p1m)/dsqrt(2d0)
            q(4,3)=(p1p+p1m)/dsqrt(2d0)

            q(1,4)=pt2x
            q(2,4)=pt2y
            q(3,4)=(p2p-p2m)/dsqrt(2d0)
            q(4,4)=(p2p+p2m)/dsqrt(2d0)

            do i=1,4
               q(i,5)=q(i,1)+q(i,2)-q(i,3)-q(i,4)
            enddo

            if(proc.eq.68)then

               elcollw=.true.
               call twobody(1,5,6,7,0d0,0d0,wt2)
               elcollw=.false.

            else

               call twojetps(mx,mq,rarr(6),rphi,ps,uh,th)

               if(proc.eq.55.or.proc.eq.62)then
                  elcollw=.true.
                  
                  if(wlp_lep)then
                  call wwmix
                  elseif(wlm_lep)then
                  call wwmix
                  endif 
                  if(wlp.eq.'mu')then
                  call twobodyw(6,8,9,0d0,mmu)
                  elseif(wlp.eq.'el')then
                  call twobodyw(6,8,9,0d0,me)
                  elseif(wlp.eq.'had')then
                  call twobodyw(6,8,9,mu_quark,md_quark)
                  else
                  call twobodyw(6,8,9,0d0,mtau)
                  endif
                  if(wlm.eq.'mu')then
                  call twobodyw(7,10,11,0d0,mmu)
                  elseif(wlm.eq.'el')then
                  call twobodyw(7,10,11,0d0,me)
                  elseif(wlm.eq.'had')then
                  call twobodyw(7,10,11,mu_quark,md_quark)
                  else
                  call twobodyw(7,10,11,0d0,mtau)
                  endif

                  elcollw=.false.
               endif

            endif

         endif

         if(diss1)call qinit(xb1,qsq1,1)
         if(diss2)call qinit(xb2,qsq2,2)

         val=wtt*wgt
         if(bin)then
            if(unw)then
            else
               call binit(val)
            endif
         endif

         if(calcmax)then
            if(wmax.lt.wtt*wgt*ren)then
               wmax=wtt*wgt*ren
               iw=iw+1
            endif
         endif

         if(unw)then
            runw=ran2()
            if(enew)then
               call unweightq(wtt*wgt*ren,runw)
            else
               call unweight(wtt*wgt*ren,runw)
            endif
         endif

 777     cs=wtt

      return
      end


