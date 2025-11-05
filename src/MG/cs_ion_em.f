ccccc EPA form factors (proton)
      subroutine cs_ionem(p,pflag,qin,out)
      implicit none
      integer i,p
      double precision pq1(4),pq2(4),xi1,xi2,xi,xii
      double precision q1(4),q2(4),qsq1,qsq2
      double precision alphaem
      double precision x1t,x2t,x1tt,x2tt,beta,xit,xbi
      double precision qsq,f1,f2,fl
      double precision out,out_tran,out_long
      double precision qi(4),pqi(4),qi2(4),qsqi2
      double precision phot,garr(-6:6),fpdf,matelem,out_t
      double precision q0,tpint,rtst,fpdf1
      REAL*8 Pmom(0:3,6)
      integer qin,pflag
      complex*16 amptest,ztt1,zout,AMP_gamq_gamgamq
      double precision f1lo,f2lo,fk1
      double precision pgg_onshell(4),pgg_diff(4),mx_onshell
      double precision px(4),pcm(4),nhel(5)
      double precision pboo(4),eta3,pmod3,cohpdf_ion,pt3
      integer i1,i2,j
      common/amptest/amptest


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
      include 'zi.f'
      include 'zoutarr.f'
      include 'gmatrices.f'
      include 'wwpars.f'
      include 'egam0.f'
      include 'pol.f'
c      include 'ewsf.f'
      include 'mion.f'
      include 'ion.f'
      include 'partonmom2.f'

      q0=0.71d0

      pol=p

cccc  qin gives initial state quark on diss side
cccc  qin = 1 : up type
cccc  qin = 2 : anti up type
cccc  qin = 3 : down type
cccc  qin = 4 : anti down type

ccccc pflag=1, photon from proton 1, pflag=2, photon from proton 2

      do i=1,4
         q1(i)=q_rf(i,1)-q_rf(i,3)
         q2(i)=q_rf(i,2)-q_rf(i,4)
c         print*,q2(i)
      enddo


      qsq1=(q_rf(4,3)-q_rf(4,1))**2-(q_rf(3,3)-q_rf(3,1))**2
     &   -(q_rf(2,3)-q_rf(2,1))**2
     &     -(q_rf(1,3)-q_rf(1,1))**2
      qsq1=-qsq1
      qsq2=(q_rf(4,4)-q_rf(4,2))**2-(q_rf(3,4)-q_rf(3,2))**2
     &   -(q_rf(2,4)-q_rf(2,2))**2
     &     -(q_rf(1,4)-q_rf(1,2))**2
      qsq2=-qsq2



c      do i=1,4
c         p1(i)=q_rf(i,1)
c         p2(i)=q_rf(i,2)
c      enddo

      xi1=-qsq1/rts/(q1(4)-q1(3))
      xi2=-qsq2/rts/(q2(4)+q2(3))

c      rtst=xi1*rts/xb1
c      xi1=xb1

c      rtst=rts

c      print*,rts,rtsnn

cccccccccccccc

      do i=1,4
         pq1(i)=0d0
         pq2(i)=0d0
      enddo

      pq1(4)=rts/2d0*xi1
      pq1(3)=pq1(4)
      pq2(4)=rts/2d0*xi2
      pq2(3)=-pq2(4)

cccccccccccccccccc

      beta=q_rf(3,1)/q_rf(4,1)

c      print*,beta,q_rf(3,2)/q_rf(4,2)
c      stop

      beta=1d0

      x1t=(q1(4)+q1(3)/beta)/rts
      x1tt=(q1(4)-q1(3)/beta)/rts
      x2t=(q2(4)-q2(3)/beta)/rts
      x2tt=(q2(4)+q2(3)/beta)/rts

c      print*,x2/x2t

ccccccccccc

      do i=1,2
c      q2(i)=0d0
      enddo
c      q2(3)=-q2(4)

c      q2(4)=dsqrt(q2(1)**2+q2(2)**2+q2(3)**2)

      do i=1,4
      pgg_onshell(i)=q1(i)+q2(i)
c      pgg_diff(i)=pgg_onshell(i)-q_rf(i,5)
c      print*,i,pgg_diff(i),pgg_onshell(i)
      enddo
      mx_onshell=dsqrt(pgg_onshell(4)**2-pgg_onshell(3)**2
     &      -pgg_onshell(2)**2-pgg_onshell(1)**2)

      do i=1,4
         px(i)=pgg_onshell(i)
         pcm(i)=p1(i)*mx_onshell/mx
      enddo

      call boost(mx_onshell,px,pcm,pboo)

c      print*,q_rf(4,6)**2-q_rf(3,6)**2-q_rf(2,6)**2-q_rf(1,6)**2
c      print*,q_rf(4,7)**2-q_rf(3,7)**2-q_rf(2,7)**2-q_rf(1,7)**2


      do i=1,4
c            print*,i,q_rf(i,6),q_rf(i,7)

c            q_rf(i,6)=pboo(i)
c            q_rf(i,7)=pgg_onshell(i)-q_rf(i,6)


c            print*,i,q_rf(i,6),q_rf(i,7)
c         q(i,6)=pboo(i)
c         q(i,7)=q(i,5)-q(i,6)
      enddo

c      print*,q_rf(4,6)**2-q_rf(3,6)**2-q_rf(2,6)**2-q_rf(1,6)**2
c      print*,q_rf(4,7)**2-q_rf(3,7)**2-q_rf(2,7)**2-q_rf(1,7)**2

c      print*,dsqrt(pgg_onshell(4)**2-pgg_onshell(3)**2
c     &-pgg_onshell(2)**2-pgg_onshell(1)**2),mx


c      do i=1,3
c            q_rf(i,6)=q_rf(i,6)+pgg_diff(i)/2d0
c            q_rf(i,7)=q_rf(i,7)+pgg_diff(i)/2d0
c      enddo
c      do i=6,7
c      q_rf(4,i)=dsqrt(q_rf(1,i)**2+q_rf(2,i)**2+q_rf(3,i)**2)
c      enddo

c       mdiff=pgg_diff(4)**2-pgg_diff(3)**2-pgg_diff(2)**2-pgg_diff(1)**2

c       xdiff=-2d0*(q_rf(4,6)*pgg_diff(4)-q_rf(3,6)*pgg_diff(3)
c      &      -q_rf(2,6)*pgg_diff(2)-q_rf(1,6)*pgg_diff(1))
c       xdiff=xdiff/mdiff

c       print*,xdiff

c       xdiff=-2d0*(q_rf(4,7)*pgg_diff(4)-q_rf(3,7)*pgg_diff(3)
c      &      -q_rf(2,7)*pgg_diff(2)-q_rf(1,7)*pgg_diff(1))
c       xdiff=xdiff/mdiff
c       xdiff=1d0-xdiff

c       print*,xdiff

cccccccccccc

c      print*,'x - ',x1t,x2t,x1,x2
c      print*,q_rf(3,1),q_rf(3,2),rts/2d0

      if(pflag.eq.2)then
         xbi=xb2
         do i=1,4
            pqi(i)=pq1(i)
            qi(i)=q1(i)
            qi2(i)=q2(i)
         enddo
         xi=xi1
         xii=xi2
         xit=x2t

         qsq=qsq1
         qsqi2=qsq2

      else
         xbi=xb1
         do i=1,4
            pqi(i)=pq2(i)
            qi(i)=q2(i)
            qi2(i)=q1(i)
         enddo
         xi=xi2
         xii=xi1
         xit=x1t

         qsq=qsq2
         qsqi2=qsq1

      endif

c      qsq=-(qi(4)**2-qi(3)**2-qi(2)**2-qi(1)**2)
c      qsqi2=-(qi2(4)**2-qi2(3)**2-qi2(2)**2-qi2(1)**2)

ccccccccccccccc

      do i=1,4
      q(i,19)=qi2(i)+pqi(i)
      enddo

      do i=1,3
         pmom(i,1)=qi2(i)
         pmom(i,2)=pqi(i)
         pmom(i,5)=pqi(i)-qi(i)
         q(i,18)=pqi(i)-qi(i)
c         q_rf(i,8)=pqi(i)-qi(i)
         pmom(i,3)=q_rf(i,6)
         pmom(i,4)=q_rf(i,7)
c         print*,i,qi2(i)
c         print*,'mom-',i,pmom(i,1),pmom(i,2),pmom(i,3)
c     &   ,pmom(i,4),pmom(i,5)
      enddo
      pmom(0,1)=qi2(4)
      pmom(0,2)=pqi(4)
      pmom(0,5)=pqi(4)-qi(4)
      q(4,18)=pqi(4)-qi(4)
c      q_rf(4,8)=pqi(4)-qi(4)
      pmom(0,3)=q_rf(4,6)
      pmom(0,4)=q_rf(4,7)
c      print*,4,qi2(4)

      call pAboost_8(1)

      pmod3=dsqrt(q(1,18)**2+q(2,18)**2+q(3,18)**2)
      eta3=(pmod3+q(3,18))/(pmod3-q(3,18))
      if(eta3.gt.0d0)then
      eta3=0.5d0*dlog(eta3)
      else
      eta3=1d10
      endif
      
      pt3=dsqrt(q(1,18)**2+q(2,18)**2)

c      if(eta3.lt.2.5d0)then
c      out=0d0
c      return
c      endif

c      if(eta3.lt.0d0)then
c      out=0d0
c      return
c      endif

c      if(eta3.lt.2.5d0)then
c      if(pt3.gt.0.5d0)then
c      out=0d0
c      return
c      endif
c      endif


c      print*,'mom-',i,pmom(0,1),pmom(0,2),pmom(0,3)
c     &   ,pmom(0,4),pmom(0,5)
c      print*,(pqi(4)-q_rf(4,6))**2-(pqi(3)-q_rf(3,6))**2
c     &      -(pqi(2)-q_rf(2,6))**2-(pqi(1)-q_rf(1,6))**2
c      print*,(pqi(4)-q_rf(4,6)-q_rf(4,7))**2
c     &      -(pqi(3)-q_rf(3,6)-q_rf(3,7))**2
c     &-(pqi(2)-q_rf(2,6)-q_rf(2,7))**2-(pqi(1)-q_rf(1,6)-q_rf(1,7))**2
c      print*,'ptgg - ',pmom(1,3)+pmom(1,4),pmom(2,3)+pmom(2,4)

      if(pmom(0,5).lt.0d0)then
      out=0d0
      return
      endif

c       if(pmom(0,5).lt.0d0)then
c       print*,-qsq1/2d0/(qi(4)-qi(3))-qi(4),pmom(0,5)
c       print*,-1d0/(qi(4)-qi(3))/2d0*
c      &(qsq1+2d0*qi(4)*(qi(4)-qi(3)))
c       print*,-1d0/(qi(4)-qi(3))/(qi(4)+qi(3))/2d0*
c      &(qsq1*(qi(4)+qi(3))+2d0*qi(4)*(qi(4)-qi(3))*(qi(4)+qi(3)))
c       print*,-1d0/(qi(4)**2-qi(3)**2)/2d0*
c      &(qsq1*(qi(4)+qi(3))+2d0*qi(4)*(qi(4)-qi(3))*(qi(4)+qi(3)))
c       print*,-1d0/(qi(4)**2-qi(3)**2)/2d0*
c      &(qsq1*(qi(4)+qi(3))+2d0*qi(4)*(qi(4)**2-qi(3)**2))
c       print*,-1d0/(qi(4)**2-qi(3)**2)/2d0*(qsq1*(qi(4)+qi(3))),
c      &      1d0/(qi(4)**2-qi(3)**2)/2d0*(2d0*qi(4)*(qi(4)**2-qi(3)**2)),
c      &      1d0/(qi(4)**2-qi(3)**2)/2d0*
c      &      (2d0*qi(4)*(qi(4)-qi(3))*(qi(4)+qi(3)))
c       print*,(qi(4)-qi(3))*(qi(4)+qi(3)),(qi(4)**2-qi(3)**2)
c       print*,''
c       print*,-1d0/(qi(4)**2-qi(3)**2)/2d0*qsq1*(qi(4)+qi(3))
c       print*,-1d0/(qi(4)-qi(3))/2d0*qsq1
c       print*,-1d0/(qi(4)**2-qi(3)**2)/2d0*2d0*qi(4)*(qi(4)**2-qi(3)**2)
c       print*,-1d0/(qi(4)-qi(3))/2d0*2d0*qi(4)*(qi(4)-qi(3))
c       print*,''
c       print*,qi(4),pqi(4),-qsq1/2d0/(qi(4)-qi(3))
c       print*,pmom(0,5)**2-pmom(1,5)**2-pmom(2,5)**2-pmom(3,5)**2
c       print*,pmom(0,5),pmom(1,5),pmom(2,5),pmom(3,5)
c       print*,pqi(4),pqi(1),pqi(2),pqi(3)
c       print*,qi(4),qi(1),qi(2),qi(3)
c       print*,qsq1,2d0*(pmom(0,5)*qi(4)-pmom(3,5)*qi(3)
c      &-pmom(2,5)*qi(2)-pmom(1,5)*qi(1))
c c      print*,qsq1,qsq2
c       print*,xi1,xb1,xi1/xb1
c       print*,xb1,qsq1/2d0/(qi(4)*q_rf(4,1)-qi(3)*q_rf(3,1))
c      &,qsq1/(qi(4)-qi(3))/rts
c       print*,qsq1*(1d0/xb1-1d0)+mp**2,mdiss1**2
c c      print*,mdiss1,mdiss2
c       print*,''
c       stop
c       endif

c$$$      qdk=(qi2(4)+beta*qi2(3))*rts/2d0
c$$$cc      do i=1,4  ! scalar pol.
c$$$         egam0(i)=q_rf(i,2)
c$$$c         +qdk/qsqi2*qi2(i)d
c$$$         egam0(i)=-egam0(i)*dsqrt(qsqi2/(qdk**2+mp**2*qsqi2))
c$$$  enddo

ccccc Correct (stable) scal pol

      egam0(4)=dsqrt(qi2(1)**2+qi2(2)**2+qi2(3)**2)/dsqrt(qsqi2)
      do i=1,3
         egam0(i)=qi2(i)*qi2(4)/dsqrt(qi2(1)**2+qi2(2)**2+qi2(3)**2)
     &        /dsqrt(qsqi2)
      enddo

      scpol=.false.

c       do i=1,5
c       nhel(i)=1
c       enddo
c       nhel(1)=-1
c       nhel(3)=-1


c       do i1=1,4
c            do i2=1,4
c            call egcalc(i1,i2)
c c           print*,i1,i2
c c           call SMATRIX_gamq_gamgamq(Pmom,MATELEM)
c            zoutarr(p,i1,i2)=AMP_gamq_gamgamq(Pmom,nhel)
c            enddo
c       enddo

c       zout=0d0
c       do i=1,1
c         zout=0d0
c         do j=1,4
c c           ztt1=zoutarr(p,i,j)*qi2(j)
c            ztt1=zoutarr(p,i,j)*q_rf(j,6)
c            if(j.lt.4)ztt1=-ztt1
c            zout=zout+ztt1
c            print*,j,zoutarr(p,i,j)
c         enddo
c         print*,i,qsqi2,zout
c c        stop
c       enddo
c       print*,''

      
      

      call SMATRIX_gamq_gamgamq(Pmom,MATELEM)

c      print*,'pmom 0 - '
c     &,pmom(0,1),pmom(0,2),pmom(0,3),pmom(0,4),pmom(0,5)
c      print*,'matelem - ',matelem

c      matelem=1d0

      

c      print*,matelem

c      if(matelem.gt.10000d0)then
c      print*,qsq,xi,xb1,xb2,matelem
c      stop
c      endif



c      matelem=1d0
         
      f2=1d0/(1d0+qsqi2/q0)**4
      f2=f2*tpint(1,dsqrt(qsqi2))**2
      f1=0d0

c      print*,tpint(1,dsqrt(qsqi2))

c      print*,'f2 - ',f2

c      f2=1d0

      fl=(1d0+4d0*mion**2*xbi**2/qsqi2)*f2-2d0*xbi*f1
      
c      fl=0d0


      out_tran=matelem

c       out_tran=out_tran*(f2*(1d0-xit/xbi+xit**2/2d0/xbi**2
c      &     +xit**2*mion**2/qsqi2)-fl*xit**2/2d0/xbi**2)   ! note divided by 2 to get averaging factor right for photons
c       out_tran=out_tran/qsqi2


      fl=(1d0+4d0*mion**2/qsqi2)*f2

      out_tran=out_tran*(f2*(1d0-xit+xit**2/2d0
     &     +xit**2*mion**2/qsqi2)-fl*xit**2/2d0)   ! note divided by 2 to get averaging factor right for photons
c      out_tran=out_tran*f2*(1d0-xit-xit**2*mion**2/qsqi2)

c      out_tran=out_tran*(q2(1)**2+q2(2)**2)/qsqi2**2*f2  ! LHL elim

c      out_tran=out_tran*(q2(1)**2+q2(2)**2)*(1d0-x2)*f2
c      out_tran=out_tran/(q2(1)**2+q2(2)**2+x2**2*mion**2)

      

c      out_tran=1d0
c      print*,xbi

ccccccc

c       scpol=.true.      

c       call SMATRIX_gamq_gamgamq(Pmom,MATELEM)

c       out_long=matelem*2d0 ! only 1 pol -> no average factor
c       out_long=out_long*(f2*(1d0-xit/xbi-xit**2*mion**2/qsqi2)
c      &     +fl*xit**2/4d0/xbi**2)

ccccccc     

      out_t=out_tran
c      out_t=out_t+out_long
      
      out_t=out_t*conv

      out_t=out_t*alphaEM(qsqi2)  ! LHL elim

      out_t=out_t/qsqi2  ! now done above

c      out_t=out_t*cohpdf_ion(xit,10d0)


c      out_t=out_t/(q2(1)**2+q2(2)**2+x2**2*mion**2)*(1d0-x2)

ccccccccccfrom rho's

      out_t=out_t/xit**2*xbi

c      print*,xbi
c      stop


cccccccc quark side PDF



      call evolvePDFphoton(xi,10d0,garr,phot)  ! TODO - scale

cccc   matrix element is for up quark so reweight accordingly

      if(qin.eq.1)then ! proton
         fpdf=(garr(2)+garr(-2)+garr(-4)+garr(-4))
         fpdf=fpdf+(garr(-1)+garr(1)+garr(-3)+garr(3))/2d0**6
         fpdf=fpdf*az

         fpdf1=(garr(1)+garr(-1)+garr(-4)+garr(-4))
         fpdf1=fpdf1+(garr(-2)+garr(2)+garr(-3)+garr(3))/2d0**6
         fpdf1=fpdf1*(an-az)

         fpdf=fpdf+fpdf1

      elseif(qin.eq.2)then ! neutron
         fpdf=(garr(1)+garr(-1)+garr(-4)+garr(-4))
         fpdf=fpdf+(garr(-2)+garr(2)+garr(-3)+garr(3))/2d0**6
         fpdf=fpdf*(an-az)
c         print*,xi
      endif

c      print*,'fpdf = ',fpdf,xi,xit

c      fpdf=1d0

      out_t=out_t*fpdf/xi

      out_t=out_t/(4d0*pi)

      out_t=out_t*4d0

cccccccccccc

      out=out_t

c      out=1d0
 
c      if((q2(1)**2+q2(2)**2).lt.0.1d0)out=0d0

c      print*,p,out,fpdf
c      print*,''
c      stop

c      out=1d0
c      print*,'out = ',out

      return
      end


