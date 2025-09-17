ccccc EPA form factors (proton)
      subroutine formfacgamoff_ion(p,q1x,q1y,q2x,q2y,zout)
      implicit none
      complex*16 zt,ztsf,zoutsf,zout
      double precision xi1,xi2,x2tt,x1tt,x1t,x2t
      double precision SF_g_ion
      double precision t1,t2,t11,t22,qsq1,qsq2
      double precision alphaem,beta
      double precision q1(2),q2(2),q1x,q2x,q1y,q2y
      integer p,i,mu,mup,nu,nup
      double precision sf1_g(4,4),sf2_g(4,4),q1p(4),q2p(4),p1(4),p2(4)
      complex*16 zw(9)
      common/zw/zw
      common/sffull/sf1_g,sf2_g

      include 'photo.f'
      include 'mp.f'
      include 'pi.f'
      include 'x.f'
      include 'xb.f'
      include 'diss.f'
      include 'vars.f'
      include 'mom.f'
      include 'ewpars.f'
      include 'proc.f'
      include 'norm.f'
      include 'ppamp.f'
      include 'zi.f'
      include 'zoutarr.f'
      include 'surv.f'
      include 'diff.f'
      include 'fbeam.f'
      include 'mion.f'
      include 'sAA.f'

c      print*,'test'

      t1=q1x**2+q1y**2
      t2=q2x**2+q2y**2

      qsq1=(q(4,3)-q(4,1))**2-(q(3,3)-q(3,1))**2-(q(2,3)-q(2,1))**2
     &     -(q(1,3)-q(1,1))**2
      qsq1=-qsq1
      qsq2=(q(4,4)-q(4,2))**2-(q(3,4)-q(3,2))**2-(q(2,4)-q(2,2))**2
     &     -(q(1,4)-q(1,2))**2
      qsq2=-qsq2

      q1(1)=q1x
      q1(2)=q1y
      q2(1)=q2x
      q2(2)=q2y


cccccccccc

      if(p.eq.1)then            ! gam-gam density matrices
         do i=1,4
            p1(i)=q(i,1)
            p2(i)=q(i,2)
         enddo

         do i=1,4
            q1p(i)=q(i,1)-q(i,3)
            q2p(i)=q(i,2)-q(i,4)
         enddo
         beta=dsqrt(1d0-4d0*mion**2/s)
         x1t=(q1p(4)+q1p(3)/beta)/rts
         x1tt=(q1p(4)-q1p(3)/beta)/rts
         x2t=(q2p(4)-q2p(3)/beta)/rts
         x2tt=(q2p(4)+q2p(3)/beta)/rts

         xi1=-qsq1/rts/(q1p(4)-q1p(3))
         xi2=-qsq2/rts/(q2p(4)+q2p(3))

c         x1t=x1
c         x2t=x2


         do mu=1,4
            do mup=1,4
               fb1=.true.
               fb2=.false.
               sf1_g(mu,mup)=SF_g_ion(diss1,mu,mup,x1t,x1tt,p1,q1p)
               fb1=.false.
               fb2=.true.
               sf2_g(mu,mup)=SF_g_ion(diss2,mu,mup,x2t,x2tt,p2,q2p)

            enddo
         enddo

      endif


cccccccc

      zoutsf=0d0

      do mu=1,4
      do mup=1,4
            do nu=1,4
            do nup=1,4

                  zt=zoutarr(p,mu,nu)*dconjg(zoutarr(p,mup,nup))

                  if(mu.lt.4)zt=-zt
                  if(mup.lt.4)zt=-zt
                  if(nu.lt.4)zt=-zt
                  if(nup.lt.4)zt=-zt

                  ztsf=zt*sf1_g(mu,mup)*sf2_g(nu,nup)
                  zoutsf=zoutsf+ztsf

               enddo
            enddo
         enddo
      enddo

c      zoutsf=zoutsf/xb1/xb2  ! multiplied by this later
      zout=dsqrt(cdabs(zoutsf))

      zout=zout*dsqrt(alphaEM(qsq1)*alphaEM(qsq2)/qsq1/qsq2)
      zout=zout*dsqrt(4d0)        ! rho normalisation

      t11=q1x**2+q1y**2
      t22=q2x**2+q2y**2

      zw(p)=zout

      return
      end
      
      
      subroutine formfacgamoff_ionp(p,q1x,q1y,q2x,q2y,zout)
      implicit none
      complex*16 zt,ztsf,zoutsf,zout
      double precision xi1,xi2,x2tt,x1tt,x1t,x2t
      double precision SF_g_ion,sf_g,SF_g_ionp
      double precision t1,t2,t11,t22,qsq1,qsq2
      double precision alphaem,beta
      double precision q1(2),q2(2),q1x,q2x,q1y,q2y
      integer p,i,mu,mup,nu,nup
      double precision sf1_g(4,4),sf2_g(4,4),q1p(4),q2p(4),p1(4),p2(4)
      complex*16 zw(9)
      common/zw/zw
      common/sffull/sf1_g,sf2_g

      include 'photo.f'
      include 'mp.f'
      include 'pi.f'
      include 'x.f'
      include 'xb.f'
      include 'diss.f'
      include 'vars.f'
      include 'mom.f'
      include 'ewpars.f'
      include 'proc.f'
      include 'norm.f'
      include 'ppamp.f'
      include 'zi.f'
      include 'zoutarr.f'
      include 'surv.f'
      include 'diff.f'
      include 'fbeam.f'
      include 'mion.f'
      include 'sAA.f'
      include 'ion_inel.f'
      include 'ion.f'

      t1=q1x**2+q1y**2
      t2=q2x**2+q2y**2

      qsq1=(q_rf(4,3)-q_rf(4,1))**2-(q_rf(3,3)-q_rf(3,1))**2
     &      -(q_rf(2,3)-q_rf(2,1))**2-(q_rf(1,3)-q_rf(1,1))**2
      qsq1=-qsq1
      qsq2=(q_rf(4,4)-q_rf(4,2))**2-(q_rf(3,4)-q_rf(3,2))**2
     &      -(q_rf(2,4)-q_rf(2,2))**2-(q_rf(1,4)-q_rf(1,2))**2
      qsq2=-qsq2

      q1(1)=q1x
      q1(2)=q1y
      q2(1)=q2x
      q2(2)=q2y


cccccccccc

      if(p.eq.1)then            ! gam-gam density matrices
         do i=1,4
            p1(i)=q_rf(i,1)
            p2(i)=q_rf(i,2)
         enddo

         do i=1,4
            q1p(i)=q_rf(i,1)-q_rf(i,3)
            q2p(i)=q_rf(i,2)-q_rf(i,4)
         enddo
         beta=dsqrt(1d0-4d0*mion**2/s)
         beta=dsqrt(1d0-4d0*mp**2/s)

         x1t=(q1p(4)+q1p(3)/beta)/rts
         x1tt=(q1p(4)-q1p(3)/beta)/rts
         x2t=(q2p(4)-q2p(3)/beta)/rts
         x2tt=(q2p(4)+q2p(3)/beta)/rts

         xi1=-qsq1/rts/(q1p(4)-q1p(3))
         xi2=-qsq2/rts/(q2p(4)+q2p(3))

         do mu=1,4
            do mup=1,4

            if(ioninel_pbeam.eq.1)then
            sf1_g(mu,mup)=SF_g_ionp(diss1,mu,mup,xb1,x1t,x1tt,p1,q1p,mx)
            sf2_g(mu,mup)=SF_g_ion(diss2,mu,mup,x2t,x2tt,p2,q2p)
            else
            sf1_g(mu,mup)=SF_g_ion(diss1,mu,mup,x1t,x1tt,p1,q1p)
            sf2_g(mu,mup)=SF_g_ionp(diss2,mu,mup,xb2,x2t,x2tt,p2,q2p,mx)
            endif

            enddo
         enddo

      endif

cccccccc

      zoutsf=0d0

      do mu=1,4
      do mup=1,4
            do nu=1,4
            do nup=1,4

                  zt=zoutarr(p,mu,nu)*dconjg(zoutarr(p,mup,nup))

                  if(mu.lt.4)zt=-zt
                  if(mup.lt.4)zt=-zt
                  if(nu.lt.4)zt=-zt
                  if(nup.lt.4)zt=-zt

                  ztsf=zt*sf1_g(mu,mup)*sf2_g(nu,nup)
                  zoutsf=zoutsf+ztsf

               enddo
            enddo
         enddo
      enddo

      zoutsf=zoutsf/xb1/xb2  ! multiplied by this later
      zout=dsqrt(cdabs(zoutsf))

      zout=zout*dsqrt(alphaEM(qsq1)*alphaEM(qsq2)/qsq1/qsq2)
      zout=zout*dsqrt(4d0)        ! rho normalisation

      t11=q1x**2+q1y**2
      t22=q2x**2+q2y**2

      if(neutron_inel)then
      zout=zout*dsqrt(an-az)
      else
      zout=zout*dsqrt(az)
      endif

      zw(p)=zout

      return
      end

      function SF_g_ionp(diss,mu,mup,xb,x,xt,p1,q1,muf)
      double precision xb,x,xt,muf,qsq1,mdiss
      double precision f1,f2,SF_g_ionp
      double precision p1(4),q1(4),q1t(4),p1t(4)
      double precision pf,qrec,crec
      integer mu,mup,i
      logical diss

      include 'gmatrices.f'
      include 'mom.f'
      include 'zi.f'
      include 'mp.f'
      include 'x.f'

      qsq1=q1(4)**2-q1(3)**2-q1(2)**2-q1(1)**2
      qsq1=-qsq1

      pf=0.25d0

      qrec=qsq1**2/4d0/mp**2+qsq1
      qrec=dsqrt(qrec)
      crec=0.75d0*qrec/pf*(1d0-(qrec/pf)**2/12d0)
      
      do i=1,2
         q1t(i)=q1(i)
      enddo
      q1t(4)=0d0
      q1t(3)=0d0

      do i=1,4
         p1t(i)=p1(i)
      enddo
      p1t(3)=-p1(3)

      mdiss=qsq1*(1d0-xb)/xb+mp**2
      mdiss=dsqrt(mdiss)
      call F1F2(diss,xb,qsq1,mdiss,f1,f2)

      if(qrec.lt.2d0*pf)then
      if(diss)then
      else
      f1=f1*crec
      f2=f2*crec
      endif
      endif

      SF_g_ionp = dble(-(d_(mu,mup))*f1+(q1t(mu)+xt*p1t(mu))*
     &     (q1t(mup)+xt*p1t(mup))*2d0*xb*f2/qsq1/x**2)

c      SF_g=SF_g*2d0                 ! In def of rho - removed as included above

      return
      end


      function SF_g_ion(diss,mu,mup,x,xt,p1,q1)
      implicit none
      double precision tpint,q0
      double precision x,xt,qsq1
      double precision f1,sf_g_ion
      double precision p1(4),q1(4),q1t(4),p1t(4)
      double precision pf,qrec,crec
      integer mu,mup,i
      logical diss

      include 'mion.f'
      include 'gmatrices.f'
      include 'mom.f'
      include 'zi.f'
      include 'mp.f'
      include 'x.f'
      include 'tppars.f'

      pf=0.25d0

      SF_g_ion=0d0
      if(mu.gt.2)return
      if(mup.gt.2)return

      qsq1=q1(4)**2-q1(3)**2-q1(2)**2-q1(1)**2
      qsq1=-qsq1

      do i=1,2
         q1t(i)=q1(i)
      enddo
      q1t(4)=0d0
      q1t(3)=0d0

      do i=1,4
         p1t(i)=p1(i)
      enddo
      p1t(3)=-p1(3)

      q0=0.71d0

      f1=1d0/(1d0+qsq1/q0)**2
      f1=f1*tpint(1,dsqrt(qsq1))

      SF_g_ion=(q1t(mu))*(q1t(mup))*4d0/qsq1/x**2
      SF_g_ion=SF_g_ion*f1**2/2d0  ! multiplied by 2 later (matches diss. norm)

      return
      end
