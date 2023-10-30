ccccc EPA form factors (proton)
      subroutine formfacgamoff(p,q1x,q1y,q2x,q2y,zout)
      implicit none
      double precision t1,t2,q1x,q2x,q1y,q2y,qsq1,qsq2
      double precision q1(2),q2(2)
      double precision xi1,xi2,x1tt,x2tt,x2t,x1t,beta
      double precision sf_g,alphaem,out,outi
      complex*16 zt,ztsf,zoutsf,zout
      integer p,i,mu,mup,nu,nup,qin
      double precision sf1_g(4,4),sf2_g(4,4),q1p(4),q2p(4),p1(4),p2(4)
      complex*16 zw(9),zoutt
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

ccccccccccc


      zoutt=0d0

      if(proc.eq.54.or.proc.eq.55)then ! for WW production only

      if(diff.eq.'dd')then  ! interpolate with full amplitude at higher Q^2
         
         if(qsq1.gt.1d0.and.qsq2.gt.1d0.and.mdiss1.gt.dsqrt(3.5d0)
     &        .and.mdiss2.gt.dsqrt(3.5d0))then

            
            call MGcross(p,out)
            
            out=out/xb1/xb2     ! normalize so consistent with above
            out=out*dsqrt(1d0-4d0*mw**2/mx**2) ! normalize so consistent with above

            zout=dsqrt(dabs(out))
            
            zout=zout*dsqrt(alphaEM(qsq1)*alphaEM(qsq2))*132.5d0

            zoutt=zout
            
            return

         elseif(qsq1.gt.1d0.and.mdiss1.gt.dsqrt(3.5d0))then
            
            out=0d0
            do qin=1,4
               call SFcalc_SD(p,2,qin,outi)
               out=out+outi
            enddo
            out=out/xb1/xb2     ! normalize so consistent with above
            out=out*dsqrt(1d0-4d0*mw**2/mx**2) ! normalize so consistent with above
            zout=dsqrt(dabs(out))

            zout=zout*dsqrt(alphaEM(qsq1)*132.5d0)
            
            return
            
         elseif(qsq2.gt.1d0.and.mdiss2.gt.dsqrt(3.5d0))then

            out=0d0
            do qin=1,4
               call SFcalc_SD(p,1,qin,outi)
               out=out+outi
            enddo
            out=out/xb1/xb2     ! normalize so consistent with above
            out=out*dsqrt(1d0-4d0*mw**2/mx**2) ! normalize so consistent with above
            zout=dsqrt(dabs(out))
            
            zout=zout*dsqrt(alphaEM(qsq2)*132.5d0) 
            
            return
            
         endif

      endif

      if(diff.eq.'sd')then

         if(diss1)then
            if(qsq1.gt.1d0.and.mdiss1.gt.dsqrt(3.5d0))then
               out=0d0
               do qin=1,4
                  call SFcalc_SD(p,2,qin,outi)
                  out=out+outi
               enddo
               out=out/xb1/xb2  ! normalize so consistent with above
               out=out*dsqrt(1d0-4d0*mw**2/mx**2) ! normalize so consistent with above
               zout=dsqrt(dabs(out))
               
               zout=zout*dsqrt(alphaEM(qsq1)*132.5d0) 

               
               return
            endif
         else
            if(qsq2.gt.1d0.and.mdiss2.gt.dsqrt(3.5d0))then
             out=0d0
               do qin=1,4
                  call SFcalc_SD(p,1,qin,outi)
                  out=out+outi
               enddo
               out=out/xb1/xb2  ! normalize so consistent with above
               out=out*dsqrt(1d0-4d0*mw**2/mx**2) ! normalize so consistent with above
               zout=dsqrt(dabs(out))

               zout=zout*dsqrt(alphaEM(qsq2)*132.5d0) 
               
               return
            endif
         endif
  
      endif

      endif
      
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
         beta=dsqrt(1d0-4d0*mp**2/s)
         x1t=(q1p(4)+q1p(3)/beta)/rts
         x1tt=(q1p(4)-q1p(3)/beta)/rts
         x2t=(q2p(4)-q2p(3)/beta)/rts
         x2tt=(q2p(4)+q2p(3)/beta)/rts

         xi1=-qsq1/rts/(q1p(4)-q1p(3))
         xi2=-qsq2/rts/(q2p(4)+q2p(3))

        
         do mu=1,4
            do mup=1,4
               fb1=.true.
               fb2=.false.
               sf1_g(mu,mup)=SF_g(diss1,mu,mup,xb1,x1t,x1tt,p1,q1p,mx)
               fb1=.false.
               fb2=.true.
               sf2_g(mu,mup)=SF_g(diss2,mu,mup,xb2,x2t,x2tt,p2,q2p,mx)         
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
      
      zw(p)=zout

      return
      end

      
      function SF_g(diss,mu,mup,xb,x,xt,p1,q1,muf)
      double precision xb,x,xt,muf,qsq1,mdiss
      double precision f1,f2,sf_g
      double precision p1(4),q1(4),q1t(4),p1t(4)
      integer mu,mup,i
      logical diss

      include 'gmatrices.f'
      include 'mom.f'
      include 'zi.f'
      include 'mp.f'
      include 'x.f'

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

      mdiss=qsq1*(1d0-xb)/xb+mp**2
      mdiss=dsqrt(mdiss)
      call F1F2(diss,xb,qsq1,mdiss,f1,f2)

      
      SF_g = -(d_(mu,mup))*f1+(q1t(mu)+xt*p1t(mu))*
     &     (q1t(mup)+xt*p1t(mup))*2d0*xb*f2/qsq1/x**2

c      SF_g=SF_g*2d0                 ! In def of rho - removed as included above                                                                  

      return
      end
