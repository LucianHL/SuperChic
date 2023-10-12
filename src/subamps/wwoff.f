ccc   gamma gamma --> l+l- subprocess amplitude - off-shell axial gauge
      subroutine wwoff_axial(p)
      implicit none
      complex*16 ampWWt,ampWWu,ampWWs,ampWWt_test,ampWW_H
      complex*16 ep(4),em(4)
      complex*16 ep_em,ep_q1,ep_q2,em_q1,em_q2
     &,pm_ep,pp_em,pp_ep,pm_em
      complex*16 ep_nL(4),em_nL(4),em2(4),ep2(4),em_em2,em_ep2,
     &     n_em,n_em2,n_ep,n_ep2,zout,em_norm,ep_norm
     &     ,ep_ep2,ztest,ampWWt_fun,ampWWu_fun,ampWWs_fun,
     &     ampaxial_fun
      double precision q1(4),q2(4),pp(4),pm(4),ap,beta,cw,sw,
     &     n_pm,n_pp,n_n,n_q1,n_q2,qsq1,time,time_old,
     &     qsq2,tl,ul,pp_q1,pp_q2,pp_pm,q1_q2,pm_q1,pm_q2
      double precision alphaEM
      integer i,j,k,l,p,mu,nu
      logical onshell
      
      include 'mom.f'
      include 'gmatrices.f'
      include 'ewpars.f'
      include 'pi.f'
      include 'vars.f'
      include 'partonmom2.f'
      include 'norm.f'
      include 'mq.f'
      include 'polww.f'
      include 'zoutarr.f'
      include 'axial.f'
      include 'wwpars.f'

      onshell=.false.
      
      beta=dsqrt(1d0-4d0*mw**2/mx**2)

      
      do i=1,4
         q1(i)=q(i,1)-q(i,3)
         q2(i)=q(i,2)-q(i,4)
         pp(i)=q(i,6)
         pm(i)=q(i,7)
      enddo

      qsq1=(q(4,3)-q(4,1))**2-(q(3,3)-q(3,1))**2-(q(2,3)-q(2,1))**2
     &     -(q(1,3)-q(1,1))**2
      qsq1=-qsq1
      qsq2=(q(4,4)-q(4,2))**2-(q(3,4)-q(3,2))**2-(q(2,4)-q(2,2))**2
     &     -(q(1,4)-q(1,2))**2
      qsq2=-qsq2

      ! LHL now done in genpol and saved
! n vector, will need to change in genpol as well if a different value is desired!
c      if (axial) then
c$$$      n(1) = 0d0
c$$$      n(2) = 0d0
c$$$      n(3) = 1d0
c$$$      n(4) = 1d0
c     endif

      if(onshell)then

         qsq1=0d0
         qsq2=0d0

         do i=1,4
            pp(i)=q(i,16)
            pm(i)=q(i,17)
            q1(i)=0d0
            q2(i)=0d0
         enddo

         q1(4)=mx/2d0
         q1(3)=mx/2d0
         q2(4)=mx/2d0
         q2(3)=-mx/2d0
         
         
      endif 
      
      do i=1,4
            
         if(p.eq.1)then         ! ++
            ep(i)=ewp(1,i)
            em(i)=ewm(1,i)
         elseif(p.eq.2)then     ! +-
            ep(i)=ewp(1,i)
            em(i)=ewm(2,i)
         elseif(p.eq.3)then     ! -+
            ep(i)=ewp(2,i)
            em(i)=ewm(1,i)
         elseif(p.eq.4)then     ! --
            ep(i)=ewp(2,i)
            em(i)=ewm(2,i)
         elseif(p.eq.5)then     ! 0+
            ep(i)=ewp(3,i)
            em(i)=ewm(1,i)
         elseif(p.eq.6)then     ! 0-
            ep(i)=ewp(3,i)
            em(i)=ewm(2,i)
         elseif(p.eq.7)then     ! +0
            ep(i)=ewp(1,i)
            em(i)=ewm(3,i)
         elseif(p.eq.8)then     ! -0
            ep(i)=ewp(2,i)
            em(i)=ewm(3,i)
         elseif(p.eq.9)then     ! 00
            ep(i)=ewp(3,i)
            em(i)=ewm(3,i)
         endif
            
      enddo


      ep_em=ep(4)*em(4)-ep(3)*em(3)-ep(2)*em(2)-ep(1)*em(1)
      ep_q1=ep(4)*q1(4)-ep(3)*q1(3)-ep(2)*q1(2)-ep(1)*q1(1)
      ep_q2=ep(4)*q2(4)-ep(3)*q2(3)-ep(2)*q2(2)-ep(1)*q2(1)
      em_q1=em(4)*q1(4)-em(3)*q1(3)-em(2)*q1(2)-em(1)*q1(1)
      em_q2=em(4)*q2(4)-em(3)*q2(3)-em(2)*q2(2)-em(1)*q2(1)
      pm_q1=pm(4)*q1(4)-pm(3)*q1(3)-pm(2)*q1(2)-pm(1)*q1(1)
      pm_q2=pm(4)*q2(4)-pm(3)*q2(3)-pm(2)*q2(2)-pm(1)*q2(1)
      pp_q1=pp(4)*q1(4)-pp(3)*q1(3)-pp(2)*q1(2)-pp(1)*q1(1)
      pp_q2=pp(4)*q2(4)-pp(3)*q2(3)-pp(2)*q2(2)-pp(1)*q2(1)
      pm_ep=pm(4)*ep(4)-pm(3)*ep(3)-pm(2)*ep(2)-pm(1)*ep(1)
      pp_em=pp(4)*em(4)-pp(3)*em(3)-pp(2)*em(2)-pp(1)*em(1)
      pp_pm=pp(4)*pm(4)-pp(3)*pm(3)-pp(2)*pm(2)-pp(1)*pm(1)
      q1_q2=q1(4)*q2(4)-q1(3)*q2(3)-q1(2)*q2(2)-q1(1)*q2(1)

      pp_ep=pp(4)*ep(4)-pp(3)*ep(3)-pp(2)*ep(2)-pp(1)*ep(1)
      pm_em=pm(4)*em(4)-pm(3)*em(3)-pm(2)*em(2)-pm(1)*em(1)

      tl=(q1(4)-pp(4))**2-(q1(3)-pp(3))**2
     &     -(q1(2)-pp(2))**2-(q1(1)-pp(1))**2
      ul=(q2(4)-pp(4))**2-(q2(3)-pp(3))**2
     &     -(q2(2)-pp(2))**2-(q2(1)-pp(1))**2

      
      n_q1 = n(4)*q1(4)-n(3)*q1(3)-n(2)*q1(2)-n(1)*q1(1)
      n_q2 = n(4)*q2(4)-n(3)*q2(3)-n(2)*q2(2)-n(1)*q2(1)
      n_pp = n(4)*pp(4)-n(3)*pp(3)-n(2)*pp(2)-n(1)*pp(1) 
      n_pm = n(4)*pm(4)-n(3)*pm(3)-n(2)*pm(2)-n(1)*pm(1)
      n_n = n(4)*n(4) - n(3)*n(3) - n(2)*n(2) - n(1)*n(1)

      do mu=1,4
         do nu=1,4

cccc  Unitary (need to change pol vectors as well though)

            if(wgauge.eq.'unitary')then
            
            ampWWt = ampWWt_fun(pp,pm,ep,em,q1,q2,ep_em,ep_q1,ep_q2,
     &           em_q1,em_q2,pp_q1,pp_q2,pm_q1,pm_q2,q1_q2,pp_pm,
     &           qsq1,qsq2,mw,d_,mu,nu,tl,ul)
            
            
            ampWWu = ampWWu_fun(pp,pm,ep,em,q1,q2,ep_em,ep_q1,ep_q2,
     &           em_q1,em_q2,pp_q1,pp_q2,pm_q1,pm_q2,q1_q2,pp_pm,
     &           qsq1,qsq2,mw,d_,mu,nu,tl,ul)
            ampWWs = ampWWs_fun(ep,em,ep_em,d_,mu,nu)
      
            
            zout=ampWWt/(tl-mw**2)+ampWWu/(ul-mw**2)+ampWWs

            else

cccc Axial
            
            zout = ampaxial_fun(pp,pm,ep,em,q1,q2,ep_em,ep_q1,ep_q2,
     &           em_q1,em_q2,pp_q1,pp_q2,pm_q1,pm_q2,tl,ul,qsq1,qsq2,
     &           mw,d_,mu,nu,pp_em,pp_ep,n,n_q1,n_q2,n_pp,n_pm,n_n)


            endif
            
            zout=zout*4d0*pi*dsqrt(alphaEM(qsq1)*alphaEM(qsq2))
            zout=zout*dsqrt(conv)
            zout=zout*dsqrt(beta)

            
cx            print*,mu,nu,zout
            
            zoutarr(p,mu,nu)=zout
            
         enddo
      enddo
            
      return
      end
