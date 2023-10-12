ccccc EPA form factors (proton)
      subroutine SFcalc_SD(p,pflag,qin,out)
      implicit none
      integer i,j,mu,mup,nup,nu,p
      double precision pq1(4),pq2(4),xi1,xi2,xi,xii
      double precision p1(4),p2(4),q1(4),q2(4),qsq1,qsq2
      double precision alphaem
      double precision x1t,x2t,x1tt,x2tt,beta,xit,xbi
      double precision sf_g,out_tot,qsq,f1,f2,fl
      double precision outn,cs_SF,out,toti,out_tran,out_long
      complex*16 ztsf,zt,zoutsf,zout_sf,sf_tran
      complex*16 zsf,asf(4,2,2),znsf,awwnsf(4,2,2)
      double precision qi(4),pqi(4),qi2(4),qsqi2,qdk
      complex*16 zamp_arr(4,2,2),zout_tot,ztot,zout_gww
      complex*16 awlh(4,2,2),zout_wq,awlhu(4,2,2),awlhd(4,2,2)
      complex*16 zout_wwg,zout_wgw
c      complex*16 asf_store(2,4,2,2),awwnsf_store(2,4,2,2)
      double precision phot,garr(-6:6),fpdf,out_lo,matelem,out_t
      double precision ap,apd,api,apdi,sw,cw,afp,afpd,eqi,eqim
      REAL*8 Pmom(0:3,6)
      integer qin,qti,iw,pflag
      logical ati
      complex*16 amptest,epg(4),emg(4)
      double precision sigtest
      double precision f1lo,f2lo,fk1
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
      include 'ewsf.f'

      
      pol=p

      sw=dsqrt(1d0-mw**2/mz**2)
      cw=dsqrt(1d0-sw**2)
      ap=1d0-aq(2)/eq(2)/sw**2
      apd=1d0-aq(1)/eq(1)/sw**2

cccc  qin gives initial state quark on diss side
cccc  qin = 1 : up type
cccc  qin = 2 : anti up type
cccc  qin = 3 : down type  
cccc  qin = 4 : anti down type

ccccc pflag=1, photon from proton 1, pflag=2, photon from proton 2      
      
      if(qin.eq.1)then
         api=ap
         apdi=apd
         ati=.false.            ! if true then antiquark
         qti=1                  ! label for keeping track of quark flavour
         eqi=eq(2)
         eqim=eq(1)
      elseif(qin.eq.2)then
         api=ap
         apdi=apd
         ati=.true.
         qti=2
         eqi=eq(2)
         eqim=eq(1)
      elseif(qin.eq.3)then
         api=apd
         apdi=ap
         ati=.false.
         qti=2
         eqi=eq(1)
         eqim=eq(2)
      else
         api=apd
         apdi=ap
         ati=.true.
         qti=1
         eqi=eq(1)
         eqim=eq(2)
      endif
      
cccccccccc

      do i=1,4
         q1(i)=q(i,1)-q(i,3)
         q2(i)=q(i,2)-q(i,4)
      enddo
         
      qsq1=(q(4,3)-q(4,1))**2-(q(3,3)-q(3,1))**2-(q(2,3)-q(2,1))**2
     &     -(q(1,3)-q(1,1))**2
      qsq1=-qsq1
      qsq2=(q(4,4)-q(4,2))**2-(q(3,4)-q(3,2))**2-(q(2,4)-q(2,2))**2
     &     -(q(1,4)-q(1,2))**2
      qsq2=-qsq2


      
      do i=1,4
         p1(i)=q(i,1)
         p2(i)=q(i,2)
      enddo

      xi1=-qsq1/rts/(q1(4)-q1(3))
      xi2=-qsq2/rts/(q2(4)+q2(3))

      
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
      
      beta=q(3,1)/q(4,1)

      x1t=(q1(4)+q1(3)/beta)/rts
      x1tt=(q1(4)-q1(3)/beta)/rts
      x2t=(q2(4)-q2(3)/beta)/rts
      x2tt=(q2(4)+q2(3)/beta)/rts

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
      endif

      qsq=-(qi(4)**2-qi(3)**2-qi(2)**2-qi(1)**2)
      qsqi2=-(qi2(4)**2-qi2(3)**2-qi2(2)**2-qi2(1)**2)

ccccccccccccccc
      
 555  do i=1,3
         pmom(i,1)=qi2(i)
         pmom(i,2)=pqi(i)
         pmom(i,3)=pqi(i)-qi(i)
         pmom(i,4)=q(i,6)
         pmom(i,5)=q(i,7)
      enddo
      pmom(0,1)=qi2(4)
      pmom(0,2)=pqi(4)
      pmom(0,3)=pqi(4)-qi(4)
      pmom(0,4)=q(4,6)
      pmom(0,5)=q(4,7)

c$$$      qdk=(qi2(4)+beta*qi2(3))*rts/2d0
c$$$cc      do i=1,4  ! scalar pol.
c$$$         egam0(i)=q(i,2)
c$$$c         +qdk/qsqi2*qi2(i)
c$$$         egam0(i)=-egam0(i)*dsqrt(qsqi2/(qdk**2+mp**2*qsqi2))
c$$$  enddo

ccccc Correct (stable) scal pol
      
      egam0(4)=dsqrt(qi2(1)**2+qi2(2)**2+qi2(3)**2)/dsqrt(qsqi2)
      do i=1,3
         egam0(i)=qi2(i)*qi2(4)/dsqrt(qi2(1)**2+qi2(2)**2+qi2(3)**2)
     &        /dsqrt(qsqi2)
      enddo

      scpol=.false.
      
      if(qin.eq.1)call SMATRIX_au(Pmom,MATELEM)
      if(qin.eq.2)call SMATRIX_aub(Pmom,MATELEM)
      if(qin.eq.3)call SMATRIX_ad(Pmom,MATELEM)
      if(qin.eq.4)call SMATRIX_adb(Pmom,MATELEM)



      if(pflag.eq.1)then
         call F1F2(diss1,xbi,qsqi2,mdiss1,f1,f2)
      else
         call F1F2(diss2,xbi,qsqi2,mdiss2,f1,f2)
      endif

c      call F1F2el(qsqi2,f1,f2)
         
      fl=(1d0+4d0*mp**2*xbi**2/qsqi2)*f2-2d0*xbi*f1

      out_tran=matelem
      out_tran=out_tran*(f2*(1d0-xit/xbi+xit**2/2d0/xbi**2
     &     +xit**2*mp**2/qsqi2)-fl*xit**2/2d0/xbi**2)   ! note divided by 2 to get averaging factor right for photons

c      out_tran=out_tran*f2*(qi2(1)**2+qi2(2)**2)/qsqi2
      
      scpol=.true.
      
      if(qin.eq.1)call SMATRIX_au(Pmom,MATELEM)
      if(qin.eq.2)call SMATRIX_aub(Pmom,MATELEM)
      if(qin.eq.3)call SMATRIX_ad(Pmom,MATELEM)
      if(qin.eq.4)call SMATRIX_adb(Pmom,MATELEM)
         
      out_long=matelem*2d0 ! only 1 pol -> no average factor
      out_long=out_long*(f2*(1d0-xit/xbi-xit**2*mp**2/qsqi2)
     &     +fl*xit**2/4d0/xbi**2)

      out_t=out_tran
      out_t=out_t+out_long
      
      
      out_t=out_t*conv
      out_t=out_t*alphaEM(qsqi2)
      out_t=out_t/qsqi2


      
      
ccccccccccfrom rho's

      out_t=out_t/xit**2*xbi
   
      
cccccccc quark side PDF

      call evolvePDFphoton(xi,dsqrt(qsq),garr,phot)

      if(qin.eq.1)then
         fpdf=(garr(2)+garr(4))
      elseif(qin.eq.2)then
         fpdf=(garr(-2)+garr(-4))
      elseif(qin.eq.3)then
         fpdf=(garr(1)+garr(3)+garr(5))
      else
         fpdf=(garr(-1)+garr(-3)+garr(-5))
      endif

      out_t=out_t*fpdf/xi

      out_t=out_t/(4d0*pi)

      out_t=out_t*4d0

ccccccccccc

      call F1F2ev(xi,qsq,f1lo,f2lo)      
      call F1F2ap(xi,qsq,f2,fl)
      f1=f2-fl
      f1=f1/2d0/xi

c      fk1=f1/f1lo
      fk1=f2/f2lo

      out_t=out_t*fk1  ! reweight so that pure SF part is ~ NNLO

cccccccccccc      
      
      out=out_t
      
      return
      end


