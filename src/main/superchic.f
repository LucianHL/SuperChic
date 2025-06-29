ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                                                         c
c     SuperChic Monte Carlo generator for central         c
c     exclusive  production.                              c
c                                                         c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     ***********************************************
c     *  v4.2                        13/03/23       *
c     *                                             *
c     *  Author: Lucian Harland-Lang                *
c     *  (l.harland-lang@ucl.ac.uk)                 *
c     *                                             *
c     *  For details see :                          *
c     *                                             *
c     *  arXiv 2506.03264 (coincident production)   *
c     *  arXiv 2303.04826 (ion dissiciation)        *
c     *  arXiv 2201.08403 (WW)                      *
c     *  arXiv 2007.12704 (v4 updates)              *
c     *  arXiv 1812.04886 (SUSY)                    *
c     *  arXiv 1810.06567                           *
c     *  arXiv 1508.02718                           *
c     *  arXiv 1405.0018 (review)                   *
c     *  arXiv 1005.0695 (quarkonia and diphoton)   *
c     *  arXiv 1011.0680 (quarkonia)                *
c     *  arXiv 1105.1626 (meson pairs)              *
c     *  arXiv 1302.2004 (meson pairs - eta/eta')   *
c     *  arXiv 1306.6661 (Skewed PDF)               *
c     *  arXiv 1306.2149 (Skewed PDF)               *
c     *                                             *
c     *  Available at :                             *
c     *  http:://projects.hepforge.org/superchic    *
c     *                                             *
c     ***********************************************
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      program superchic
      implicit none
      double precision randum,ran2,t2
      double precision sd,sdo,sd1,chi2a,avgio,avgi1,avgi
      double precision mass,mup,ms,md,beta
      integer i,j,k,itmx1
      integer nhistmax
      integer outl
      integer iinc,ncallu
      logical histol
      character*100 dum
      integer idum
      COMMON /ranno/ idum

      include 'pdfinf.f'
      include 'genunw.f'
      include 'survin.f'
      include 'mesflag.f'
      include 'procn.f'
      include 'gencuts.f'
      include 'pi.f'
      include 'pdg.f'
      include 'unweighted.f'
      include 'surv.f'
      include 'bin.f'
      include 'ewpars.f'
      include 'zi.f'
      include 'vars.f'
      include 'survpars.f'
      include 'mom.f'
      include 'range.f'
      include 'polarization.f'
      include 'proc.f'
      include 'mq.f'
      include 'norm.f'
      include 'mixing.f'
      include 'meta.f'
      include 'mres.f'
      include 'forward.f'
      include 'mfact.f'
      include 'jetalg.f'
      include 'decay.f'
      include 'record.f'
      include 'mp.f'
      include 'quarkonia.f'
      include 'scorr.f'
      include 'output.f'
      include 'nhist.f'
      include 'intag.f'
      include 'vegas.f'
      include 'vegaspars.f'
      include 'hepevt.f'
      include 'leshouches.f'
      include 'photo.f'
      include 'nsurv.f'
      include 'eff.f'
      include 'prec.f'
      include 'wmax.f'
      include 'wtmax.f'
      include 'mkp.f'
      include 'mpip.f'
      include 'widths.f'
      include 'lep.f'
      include 'monopar.f'
      include 'gax.f'
      include 'beam.f'
      include 'ion.f'
      include 'mion.f'
      include 'ionqcd.f'
      include 'qcd.f'
      include 'mn.f'
      include 'rech.f'
      include 'ptXcuts.f'
      include 'mcharg.f'
      include 'varsi.f'
      include 'spA.f'
      include 'sAA.f'
      include 'diff.f'
      include 'diss.f'
      include 'enew.f'
      include 'gamma.f'
      include 'wwpars.f'
      include 'wdecay.f'
      include 'p0Xn.f'
      include 'mxs.f'
      include 'tau.f'
      character*10 tdiff,tbeam,bp

      call EXECUTE_COMMAND_LINE('mkdir -p inputs evrecs outputs')

ccccccc

      read(*,*)dum
      read(*,*)dum
      read(*,*)dum
      read(*,*)rts
      read(*,*)isurv
      read(*,*)intag
      read(*,*)dum
      read(*,*)dum
      read(*,*)PDFname
      read(*,*)PDFmember
      SFPDFname='SF_MSHT20qed_nnlo'
      SFPDFmember=0
      read(*,*)dum
      read(*,*)proc
      read(*,*)beam
      read(*,*)outtag
      read(*,*)sfaci
      read(*,*)dum
      read(*,*)dum
      read(*,*)dum
      read(*,*)diff
c      read(*,*)elcoll
      read(*,*)dum
      read(*,*)an
      read(*,*)az
      read(*,*)rz
      read(*,*)dz
      read(*,*)rn
      read(*,*)dn
      read(*,*)ionqcd
      read(*,*)ionbreakup
      read(*,*)fAA
      read(*,*)fracsigX
      read(*,*)wrho
      read(*,*)yrho
      read(*,*)accrho
      read(*,*)dum
      read(*,*)dum
      read(*,*)dum
      read(*,*)ncall
      read(*,*)itmx
      read(*,*)prec
      read(*,*)ncall1
      read(*,*)inccall
      read(*,*)itend
      read(*,*)iseed
      read(*,*)dum
      read(*,*)dum
      read(*,*)dum
      read(*,*)genunw
      read(*,*)nev
      read(*,*)erec
      read(*,*)readwt
      read(*,*)wtmax
      read(*,*)dum
      read(*,*)dum
      read(*,*)dum
      read(*,*)ymin
      read(*,*)ymax
      read(*,*)mmin
      read(*,*)mmax
      read(*,*)gencuts
      read(*,*)scorr
      read(*,*)fwidth
      read(*,*)dum
      read(*,*)dum
      read(*,*)dum
      read(*,*)dum
      read(*,*)dum
      read(*,*)ptxmax
      read(*,*)dum
      read(*,*)dum
      read(*,*)dum
      read(*,*)ptamin
      read(*,*)ptbmin
      read(*,*)etaamin
      read(*,*)etaamax
      read(*,*)etabmin
      read(*,*)etabmax
      read(*,*)acoabmax
      read(*,*)dum
      read(*,*)dum
      read(*,*)dum
      read(*,*)ptamin3
      read(*,*)ptbmin3
      read(*,*)ptcmin3
      read(*,*)etaamin3
      read(*,*)etaamax3
      read(*,*)etabmin3
      read(*,*)etabmax3
      read(*,*)etacmin3
      read(*,*)etacmax3
      read(*,*)dum
      read(*,*)dum
      read(*,*)dum
      read(*,*)ptamin4
      read(*,*)ptbmin4
      read(*,*)ptcmin4
      read(*,*)ptdmin4
      read(*,*)etaamin4
      read(*,*)etaamax4
      read(*,*)etabmin4
      read(*,*)etabmax4
      read(*,*)etacmin4
      read(*,*)etacmax4
      read(*,*)etadmin4
      read(*,*)etadmax4
      read(*,*)dum
      read(*,*)dum
      read(*,*)dum
      read(*,*)ptamin6
      read(*,*)ptbmin6
      read(*,*)ptcmin6
      read(*,*)ptdmin6
      read(*,*)ptemin6
      read(*,*)ptfmin6
      read(*,*)etaamin6
      read(*,*)etaamax6
      read(*,*)etabmin6
      read(*,*)etabmax6
      read(*,*)etacmin6
      read(*,*)etacmax6
      read(*,*)etadmin6
      read(*,*)etadmax6
      read(*,*)etaemin6
      read(*,*)etaemax6
      read(*,*)etafmin6
      read(*,*)etafmax6
      read(*,*)dum
      read(*,*)dum
      read(*,*)dum
      read(*,*)rjet
      read(*,*)jalg
      read(*,*)dum
      read(*,*)dum
      read(*,*)dum
      read(*,*)m2b
      read(*,*)pdgid(6)
      read(*,*)pdgid(7)
      read(*,*)dum
      read(*,*)dum
      read(*,*)dum
      read(*,*)malp
      read(*,*)gax
      read(*,*)alpt
      read(*,*)dum
      read(*,*)dum
      read(*,*)dum
      read(*,*)mpol
      read(*,*)mmon
      read(*,*)gamm
      read(*,*)dum
      read(*,*)dum
      read(*,*)dum
      read(*,*)mcharg
      read(*,*)mneut
      read(*,*)dum
      read(*,*)dum
      read(*,*)dum
      read(*,*)wlp
      read(*,*)wlm
      read(*,*)dum
      read(*,*)dum
      read(*,*)dum
      read(*,*)tau
      read(*,*)mxs
      read(*,*)dum
      read(*,*)dum
      read(*,*)dum
      read(*,*)atau
      read(*,*)dtau
      read(*,*)
      read(*,*)calc_tau_coeff
      read(*,*)tau_mom
      read(*,*)tau_coeff

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      if(beam.ne.'ion')wrho=.false.
      if(wrho)ionbreakup=.true.
      if(fAA.eq.'10')fAA='01'
      if(fAA.eq.'1A')fAA='A1'
      if(fAA.eq.'XA')fAA='AX'
      if(fAA.eq.'0A')fAA='A0'
      if(fAA.eq.'X0')fAA='0X'
      if(fAA.eq.'X1')fAA='1X'
      
      int_01=.false.
      if(fAA.eq.'01'.or.fAA.eq.'A1'.or.fAA.eq.'AX')int_01=.true.
      if(wrho)then
         if(fAA.eq.'AA')int_01=.true.
         if(fAA.eq.'01')int_01=.false.
         if(fAA.eq.'AX')int_01=.false.
         if(fAA.eq.'A1')int_01=.false.
         if(fAA.eq.'A0')int_01=.true.
         if(fAA.eq.'00')int_01=.true.
      endif

      tdiff=diff
      if(diff.eq.'sda'.or.diff.eq.'sdb')then
         tdiff='sd'
      endif
      tbeam=beam
      if (beam.eq.'el') tbeam='ee'
      if (beam.eq.'ionp') tbeam='pA'
      if (beam.eq.'ion') tbeam='AA'
      if (beam.eq.'prot') tbeam='pp'
      bp=tbeam(1:len(trim(tbeam)))//tdiff(1:len(trim(tdiff)))

      if(proc.eq.1.and.(bp.eq.'ppel')) goto 111
      if(proc.eq.2.and.(bp.eq.'ppel')) goto 111
      if(proc.eq.3.and.(bp.eq.'ppel')) goto 111
      if(proc.eq.4.and.(bp.eq.'ppel')) goto 111
      if(proc.eq.5.and.(bp.eq.'ppel')) goto 111
      if(proc.eq.6.and.(bp.eq.'ppel')) goto 111
      if(proc.eq.7.and.(bp.eq.'ppel')) goto 111
      if(proc.eq.8.and.(bp.eq.'ppel')) goto 111
      if(proc.eq.9.and.(bp.eq.'ppel')) goto 111
      if(proc.eq.10.and.(bp.eq.'ppel')) goto 111
      if(proc.eq.11.and.(bp.eq.'ppel')) goto 111
      if(proc.eq.12.and.(bp.eq.'ppel')) goto 111
      if(proc.eq.13.and.(bp.eq.'ppel')) goto 111
      if(proc.eq.14.and.(bp.eq.'ppel')) goto 111
      if(proc.eq.15.and.(bp.eq.'ppel')) goto 111
      if(proc.eq.16.and.(bp.eq.'ppel')) goto 111
      if(proc.eq.17.and.(bp.eq.'ppel')) goto 111
      if(proc.eq.18.and.(bp.eq.'ppel')) goto 111
      if(proc.eq.19.and.(bp.eq.'ppel')) goto 111
      if(proc.eq.20.and.(bp.eq.'ppel')) goto 111
      if(proc.eq.21.and.(bp.eq.'ppel')) goto 111
      if(proc.eq.22.and.(bp.eq.'ppel')) goto 111
      if(proc.eq.23.and.(bp.eq.'ppel')) goto 111
      if(proc.eq.24.and.(bp.eq.'ppel')) goto 111
      if(proc.eq.25.and.(bp.eq.'ppel')) goto 111
      if(proc.eq.26.and.(bp.eq.'ppel')) goto 111
      if(proc.eq.27.and.(bp.eq.'ppel')) goto 111
      if(proc.eq.28.and.(bp.eq.'ppel')) goto 111
      if(proc.eq.29.and.(bp.eq.'ppel')) goto 111
      if(proc.eq.30.and.(bp.eq.'ppel')) goto 111
      if(proc.eq.31.and.(bp.eq.'ppel')) goto 111
      if(proc.eq.32.and.(bp.eq.'ppel')) goto 111
      if(proc.eq.33.and.(bp.eq.'ppel')) goto 111
      if(proc.eq.34.and.(bp.eq.'ppel')) goto 111
      if(proc.eq.35.and.(bp.eq.'ppel')) goto 111
      if(proc.eq.36.and.(bp.eq.'ppel')) goto 111
      if(proc.eq.37.and.(bp.eq.'ppel')) goto 111
      if(proc.eq.38.and.(bp.eq.'ppel')) goto 111
      if(proc.eq.39.and.(bp.eq.'ppel')) goto 111
      if(proc.eq.40.and.(bp.eq.'ppel')) goto 111
      if(proc.eq.41.and.(bp.eq.'ppel')) goto 111
      if(proc.eq.42.and.(bp.eq.'ppel')) goto 111
      if(proc.eq.43.and.(bp.eq.'ppel')) goto 111
      if(proc.eq.44.and.(bp.eq.'ppel')) goto 111
      if(proc.eq.45.and.(bp.eq.'ppel')) goto 111
      if(proc.eq.46.and.(bp.eq.'ppel')) goto 111
      if(proc.eq.47.and.(bp.eq.'ppel')) goto 111
      if(proc.eq.48.and.(bp.eq.'pAel'.or.bp.eq.'ppel')) goto 111
      if(proc.eq.49.and.(bp.eq.'pAel'.or.bp.eq.'ppel')) goto 111
      if(proc.eq.50.and.(bp.eq.'pAel'.or.bp.eq.'ppel')) goto 111
      if(proc.eq.51.and.(bp.eq.'pAel'.or.bp.eq.'ppel')) goto 111
      if(proc.eq.52.and.(bp.eq.'pAel'.or.bp.eq.'ppel')) goto 111
      if(proc.eq.53.and.(bp.eq.'pAel'.or.bp.eq.'ppel')) goto 111
      if(proc.eq.54.and.(bp.eq.'pAel')) goto 111
      if(proc.eq.55.and.(bp.eq.'pAel'.or.bp.eq.'AAel'.or.bp.eq.'eeel'
     & .or.bp.eq.'ppdd'.or.bp.eq.'ppsd'.or.bp.eq.'ppel')) goto 111
      if(proc.eq.56.and.(bp.eq.'pAel'.or.bp.eq.'AAel'.or.bp.eq.'eeel'
     &.or.bp.eq.'ppdd'.or.bp.eq.'ppsd'.or.bp.eq.'ppel')) goto 111
      if(proc.eq.57.and.(bp.eq.'pAel'.or.bp.eq.'AAel'.or.bp.eq.'eeel'
     &.or.bp.eq.'ppdd'.or.bp.eq.'ppsd'.or.bp.eq.'ppel')) goto 111
      if(proc.eq.58.and.(bp.eq.'pAel'.or.bp.eq.'AAel'.or.bp.eq.'eeel'
     &.or.bp.eq.'ppdd'.or.bp.eq.'ppsd'.or.bp.eq.'ppel')) goto 111
      if(proc.eq.59.and.(bp.eq.'pAel'.or.bp.eq.'AAel'.or.bp.eq.'eeel'
     &.or.bp.eq.'ppel')) goto 111
      if(proc.eq.60.and.(bp.eq.'pAel'.or.bp.eq.'AAel'.or.bp.eq.'eeel'
     &.or.bp.eq.'ppel')) goto 111
      if(proc.eq.61.and.(bp.eq.'pAel')) goto 111
      if(proc.eq.68.and.(bp.eq.'pAel'.or.bp.eq.'AAel'.or.bp.eq.'eeel'
     &.or.bp.eq.'ppdd'.or.bp.eq.'ppsd'.or.bp.eq.'ppel')) goto 111
      if(proc.eq.69.and.(bp.eq.'pAel'.or.bp.eq.'AAel'.or.bp.eq.'eeel'
     &.or.bp.eq.'ppel')) goto 111
      if(proc.eq.70.and.(bp.eq.'pAel'.or.bp.eq.'AAel'.or.bp.eq.'eeel'
     & .or.bp.eq.'ppel')) goto 111
      if(proc.eq.71.and.(bp.eq.'pAel'.or.bp.eq.'AAel'.or.bp.eq.'eeel'
     &.or.bp.eq.'ppel')) goto 111
      if(proc.eq.72.and.(bp.eq.'pAel'.or.bp.eq.'AAel'.or.bp.eq.'eeel'
     & .or.bp.eq.'ppel')) goto 111
      if(proc.eq.73.and.(bp.eq.'pAel'.or.bp.eq.'AAel'.or.bp.eq.'eeel'
     &.or.bp.eq.'ppel')) goto 111
      if(proc.eq.74.and.(bp.eq.'pAel'.or.bp.eq.'AAel'.or.bp.eq.'eeel'
     & .or.bp.eq.'ppel')) goto 111
      if(proc.eq.75.and.(bp.eq.'pAel'.or.bp.eq.'AAel'.or.bp.eq.'eeel'
     &.or.bp.eq.'ppel')) goto 111
      if(proc.eq.76.and.(bp.eq.'pAel'.or.bp.eq.'AAel'.or.bp.eq.'eeel'
     &.or.bp.eq.'ppel')) goto 111
      if(proc.eq.77.and.(bp.eq.'pAel'.or.bp.eq.'ppel')) goto 111
      if(proc.eq.78.and.(bp.eq.'pAel'.or.bp.eq.'ppel')) goto 111
      if(proc.eq.79.and.(bp.eq.'pAel'.or.bp.eq.'ppel')) goto 111
      if(proc.eq.82.and.(bp.eq.'pAel'.or.bp.eq.'eeel'
     &.or.bp.eq.'ppel')) goto 111
      if(proc.eq.83.and.(bp.eq.'pAel'.or.bp.eq.'eeel'
     &.or.bp.eq.'ppel')) goto 111
      if(proc.eq.84.and.(bp.eq.'pAel'.or.bp.eq.'eeel'
     &.or.bp.eq.'ppel')) goto 111
      if(proc.le.47.and.beam.eq.'ion')goto 111

      write(*,*)'Usupported process and beam combination'
      write(*,*)'bp->',bp,'<-'
      write(*,*)'proc->',proc,'<-'
      STOP 1
 111  continue 
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      approx=.false.

ccccccccccccc

      wpol='tot'
      ftype='tot'
      higgs=.true.
      sfonly=.false.
      addnsf=.true.
      wgauge='axial'

ccccccccc

      if(gencuts.eqv..false.)then
         print*,'WARNING - gencuts=.false., no cuts on decay products'
      endif

      if(sfaci.eqv..false.)then
         if(ionbreakup)then
            print*,'sfaci=.false. -> ionbreakup set to .false'
         endif
      endif

      if(sfaci)then
         if(proc.eq.82.or.proc.eq.83.or.proc.eq.84)then
            print*,"Simplified model -> setting sfaci to false"
            sfaci=.false.
         endif
      endif

      diffsd=diff
      if(diff.eq.'sda'.or.diff.eq.'sdb')then
         diff='sd'
      else
         diffsd='n'
      endif

      offshell=.false.

      forward=.false.

      call length(outtag,outl)


      open(45,file='evrecs/evrec'//outtag(1:outl)//'.dat')
      wmax=0d0
      evnum=0

      if(genunw)then
      else
         readwt=.false.
      endif
      if(readwt)wmax=wtmax

      if(erec.eq.'hepmc')then
         erech=.true.
         erec='lhe'
      endif

      iw=0

      gf=1.16639d-5
      v=dsqrt(1d0/dsqrt(2d0)/gf)
      mt=173d0
      mb=4.75d0
      mc=1.4d0
      mmu=0.10566d0
      mpsi=3.096916d0
      mpsip=3.686109d0
      mups=9.46030d0
      mchic0=3.41475d0
      mchib0=9.85944d0
      mp=0.938272046d0
      mn=0.939565413d0
c      mwx=80.318d0
      mw=80.419d0
      me=0.511d-3
      mtau=1.77682d0
      mpip=0.13957018d0
      mkp=0.493677d0
      alpha=7.2974d-3

      pi=dacos(-1d0)
      conv=389379d3
      zi=(0d0,1d0)
      mup=0.062d0
      md=0.083d0
      ms=0.215d0
      md_quark=md
      mu_quark=mup


      NPlin=.false.  ! if true only linear terms in a_tau,d_tau included (not available by default)
      if(NPlin)int_atauonly=.false.

      if(proc.ne.58)then
         print*
     &,'[calc_tau_coeff] can only be true for tau tau - set to false'
         calc_tau_coeff=.false.
      endif

      if(calc_tau_coeff)then
         if(tau_mom.eq.'atau')then
            atau=1d0
            dtau=0d0
         elseif(tau_mom.eq.'dtau')then
            atau=0d0
            dtau=1d-15
         else  
            print*,'Incorrect calc_tau_coeff flag -- STOP'
            stop
         endif 
         if(tau_coeff.eq.0)then
            atau=0d0
            dtau=0d0
            deltau=.false.
            atau_only=.false.
            atau_lin=.false.
            atau_quad=.false.
            int_atauonly=.false.
         elseif(tau_coeff.eq.1)then
            deltau=.true.
            atau_only=.false.
            atau_lin=.true.
            atau_quad=.false.
            int_atauonly=.false.
         elseif(tau_coeff.eq.2)then
            deltau=.true.
            atau_only=.false.
            atau_lin=.false.
            atau_quad=.true.
            int_atauonly=.false.
         elseif(tau_coeff.eq.3)then
            deltau=.false.
            atau_only=.true.
            atau_lin=.false.
            atau_quad=.false.
            int_atauonly=.true.
         elseif(tau_coeff.eq.3)then
            deltau=.false.
            atau_only=.true.
            atau_lin=.false.
            atau_quad=.false.
            int_atauonly=.true.
         elseif(tau_coeff.eq.4)then
            deltau=.false.
            atau_only=.true.
            atau_lin=.false.
            atau_quad=.true.
            int_atauonly=.false.
         else  
            print*,'Incorrect tau_coeff flag -- STOP'
            stop
         endif 
         if(tau_mom.eq.'dtau')then
            if(tau_coeff.eq.1.or.tau_coeff.eq.3)then
               print*,'Only even tau_coeff non-zero for dtau -- STOP'
               stop
            endif
         endif
      else  
         deltau=.false.
         atau_only=.false.
         atau_lin=.false.
         atau_quad=.false.
         int_atauonly=.false.
      endif 


      if(deltau)atau_only=.false.
      dtau=dtau/0.1973d-13*2d0*dsqrt(pi/1.325070D+02) ! convert to GeV^-1
      GC_6515=dsqrt(pi/1.325070D+02)/mtau*atau/2d0*zi
      GC_6506=-dtau/2d0

      rmf1( 1) = 1d-10
      rmf1( 2) = me
      rmf1( 3) = 1d-10
      rmf1( 4) = mmu
      rmf1( 5) = 1d-10
      rmf1( 6) = mtau
      rmf1( 7) = 0.062d0
      rmf1( 8) = 0.083d0
      rmf1( 9) = mc
      rmf1(10) = 0.215d0
      rmf1(11) = mt
      rmf1(12) = mb

      rmf1( 1) = me
      rmf1( 2) = mmu
      rmf1( 3) = mtau
      rmf1( 4) = md
      rmf1( 5) = mup
      rmf1( 6) = ms
      rmf1( 7) = mc
      rmf1( 8) = mb
      rmf1( 9) = mt


      mq=0d0
      hel=1
      mes=.false.
      mfact='mx'
      forward=.false.
      decay2=.false.
      decay3=.false.
      decay4=.false.
      decay6=.false.

cccccccccccc

      do i=1,20
         jdahep(1,i)=0
         jdahep(2,i)=0
      enddo

cccccccccccccccccccccccccc

      call inpdf
      call supinit

      if(proc.eq.54.or.proc.eq.55.or.proc.eq.58)then
         call setpara('param_card.dat') !set parameters for MG calculation
      endif

      s2int=8
      if(beam.eq.'ionp')s2int=16
      if(diff.eq.'el'.and.gamma.eqv..true.)s2int=16


      call header
      call gaminit
      call gaminit_comb
      call readcoh
      call gaussinit

cccccccccccccccccccccccccc

      if(diff.eq.'sd'.or.diff.eq.'dd')then
         if(offshell.eqv..false.)then
            print*,'Dissociation not currently supported'//
     &           ' for this process/beam - STOP'
            STOP 1
         endif
         if(erec.eq.'hepevt'.or.erec.eq.'hepmc')then
            print*,'Dissociation currently only supported with LHE'
         endif
      endif

      elcoll=.false.
      difftot=.false.
      if(diff.eq.'el')then
         diss1=.false.
         diss2=.false.
      elseif(diff.eq.'sd')then
         ndim=ndim+1
         if(genunw)elcoll=.true.
      elseif(diff.eq.'dd')then
         diss1=.true.
         diss2=.true.
         ndim=ndim+2
      elseif(diff.eq.'tot')then
         ndim=ndim+3
         difftot=.true.
      endif

      if(diff.eq.'tot'.or.diff.eq.'sd'.or.diff.eq.'dd')then
         call apfelinit
         call calcs2diss
      endif

cccccccccccccccccccccccccc

      if(mes)then
         call calcmes
         call wfinit
      endif

ccccccccccccccccccccccccc

      if(beam.eq.'el')then
         if(sfaci)then
            print*,'Error : must have sfaci = .false. for initial-state
     &electrons'
            STOP 1
         endif
      endif

ccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccc

      s=rts**2

      if(beam.eq.'prot')then
         beta=dsqrt(1d0-4d0*mp**2/s)
      elseif(beam.eq.'el')then
         beta=dsqrt(1d0-4d0*me**2/s)
      elseif(beam.eq.'ionp'.or.beam.eq.'ion')then
c     mion=mp*an
         mion=mp*az+(an-az)*mn
         rtsi=rts
         si=s
      endif

      q(1,1)=0d0
      q(2,1)=0d0
      q(3,1)=rts/2d0*beta
      q(4,1)=rts/2d0

      q(1,2)=0d0
      q(2,2)=0d0
      q(3,2)=-rts/2d0*beta
      q(4,2)=rts/2d0

      if(beam.eq.'ionp')call pAinit
      if(beam.eq.'ion')call AAinit

      if(beam.eq.'prot')then
         pdgid(1)=2212
         pdgid(2)=2212
         pdgid(3)=2212
         pdgid(4)=2212
      elseif(beam.eq.'el')then
         pdgid(1)=11
         pdgid(2)=-11
         pdgid(3)=11
         pdgid(4)=-11
      elseif(beam.eq.'ion')then
         pdgid(1)=1000000000
         pdgid(1)=pdgid(1)+nint(az)*10000
         pdgid(1)=pdgid(1)+nint(an)*10
         pdgid(2)=pdgid(1)
         pdgid(3)=pdgid(1)
         pdgid(4)=pdgid(1)
      elseif(beam.eq.'ionp')then
         pdgid(2)=1000000000
         pdgid(2)=pdgid(2)+nint(az)*10000
         pdgid(2)=pdgid(2)+nint(an)*10
         pdgid(1)=2212
         pdgid(3)=pdgid(1)
         pdgid(4)=pdgid(2)
      endif

ccccccccccccccccccccccccccccccccccccccccccccccc
cccc     HEPEVT
ccccccccccccccccccccccccccccccccccccccccccccccc

      do k=1,2
         do j=1,4
            phep(j,k)=q(j,k)
         enddo
         phep(5,k)=mass(k)
         if(beam.eq.'el')then
            phep(5,k)=me
         elseif(beam.eq.'prot')then
            phep(5,k)=mp
         elseif(beam.eq.'ion'.or.beam.eq.'ionp')then
            phep(5,k)=mion
         endif
      enddo
      do k=1,20
         do j=1,4
            vhep(j,k)=0d0
         enddo
      enddo
      isthep(1)=4
      isthep(2)=4
      isthep(3)=1
      isthep(4)=1
      jmohep(1,5)=1
      jmohep(2,5)=2

      jdahep(1,1)=3
      jdahep(2,1)=5
      jdahep(1,2)=3
      jdahep(2,2)=5
      jdahep(1,3)=0
      jdahep(2,3)=0
      jdahep(1,4)=0
      jdahep(2,4)=0


ccccccccccccccccccccccccccccccccccccccccccccccc
ccccc Les Houches
ccccccccccccccccccccccccccccccccccccccccccccccc

      nprup=1
      idwtup=3
      pdfsup(1)=0
      pdfsup(2)=0
      pdfgup(1)=0
      pdfgup(2)=0
      idprup=0
      xwgtup=1d0
      aqedup=alpha

ccc   NEW LHE init

      pdfsup(1)=247000
      pdfsup(2)=247000
      pdfgup(1)=0
      pdfgup(2)=0

      do k=1,2
         do j=1,4
            pup(j,k)=q(j,k)
         enddo
         pup(5,k)=mass(k)
         if(beam.eq.'el')then
            pup(5,k)=me
         elseif(beam.eq.'prot')then
            pup(5,k)=mp
         elseif(beam.eq.'ion'.or.beam.eq.'ionp')then
            pup(5,k)=mion
         endif
      enddo
      istup(1)=-1
      istup(2)=-1
      istup(3)=1
      istup(4)=1
      mothup(1,1)=0
      mothup(2,1)=0
      mothup(1,2)=0
      mothup(2,2)=0

      icolup(1,1)=0
      icolup(2,1)=0
      icolup(1,2)=0
      icolup(2,2)=0
      icolup(1,3)=0
      icolup(2,3)=0
      icolup(1,4)=0
      icolup(2,4)=0
      icolup(1,5)=0
      icolup(2,5)=0
      do i=1,20
         vtimup(i)=0
         spinup(i)=9
      enddo


ccccccccc

      if(proc.eq.18.or.proc.eq.19.or.proc.eq.20)then
         nup=11
      elseif(proc.eq.83.or.proc.eq.84)then
         nup=9
      elseif(proc.eq.54.or.proc.eq.55)then
         nup=11
      elseif(proc.eq.73.or.proc.eq.74.or.proc.eq.75)then
         nup=13
      elseif(proc.eq.76)then
         nup=11
      elseif(decay4)then
         if(proc.eq.53)then
            nup=10
         else
            nup=9
         endif
      elseif(decay6)then
         nup=11
      elseif(dps.eq.2.or.decay2.or.dps.eq.12)then
         nup=7
      elseif(dps.eq.3)then
         nup=8
      elseif(decay3)then
         nup=9
      elseif(dps.eq.1)then
         nup=5
      endif


      if(diff.eq.'sd')nup=nup-1
      if(diff.eq.'el')nup=nup-2

ccccccccc

      surv=1d0
      if( (beam.eq.'prot') .or. 
     & (beam.eq.'ion' .and. (ionqcd.eq.'coh'.or.ionqcd.eq.'incoh')) .or.
     & (beam.eq.'ionp' .and. 
     & (ionqcd.eq.'coh'.or.ionqcd.eq.'incoh'))) then
         call initparsr(isurv)
         call readscreen
         if(beam.eq.'prot'.or.ionqcd.eq.'incoh')surv=1d0/norm**2
      endif

      if(qcd)then
C         call calcsud
C         call calchg
         call readsud
         call readhg
      endif



      if(beam.eq.'ion'.or.beam.eq.'ionp')call ioninit


      if(beam.eq.'ionp')then
         rts=rtspa
         s=spa
      elseif(beam.eq.'ion')then
         rts=rtsaa
         s=saa
      endif

cccccccccccc

      nhist=0
      nhistmax=20

ccccccccccccccc

      histol=.true.

ccccccc    initialise histograms

      if(histol)call inithist(nhistmax)

cccccccccccccccc

      neff=0
      neff0=0

      do i=1,10
         xu(i)=1d0
         xl(i)=0d0
      enddo

      ACC=-1D0
      NPRN=1

      idum=-abs(iseed)
      randum=ran2()


      ITMX1=1

      bin=.false.
      sfac=.false.
      unw=.false.
      calcmax=.false.
      iinc=1


      print*,''
      print*,'**********************************************************
     &**************'
      print*,'                Vegas: initialisation run                '
      print*,'                Note : outputs *bare* cross section only '
      print*,'**********************************************************
     &**************'

      CALL VEGAS(AVGI,SD,CHI2A)
c      CALL VEGAS(cs,AVGI,SD,CHI2A)

      if(readwt)goto 779

      print*,''
      print*,'**********************************************************
     &**************'
      print*,'                Vegas : main run                '
      print*,'**********************************************************
     &**************'

      if(gencuts)then

      print*,''
      print*,'**********************************************************
     &**************'
      print*,''
      write(6,19)dble(neff)/dble(neff0)*100d0
      print*,''

      iinc=nint(dble(neff0)/dble(neff))

      write(6,20)iinc
      print*,''
      print*,'**********************************************************
     &**************'
      print*,''

 19   FORMAT(' Cut Efficiency = ',G10.4,'%')
 20   FORMAT(' => multiply NCALL by ',i4)

      endif


      ITMX=ITMX1
      NCALL=NCALL1
      avgi1=avgi
      sd1=sd

      ncall=ncall*iinc
      inccall=inccall*iinc


 779  bin=.true.

      ncallu=ncall

      sfac=sfaci
      unw=.false.
      calcmax=.true.
      ren=1d0

      it=1
      itmx=1



      if(readwt)goto 778

      ren=dble(ncall)

      CALL VEGAS1(AVGI,SD,CHI2A)
c      CALL VEGAS1(cs,AVGI,SD,CHI2A)

      prec=prec*1d-2

c      open(40,file="phitesta.dat")

c      write(40,*)7
c      print*,'test'


 777  if(dabs(sd).gt.dabs(avgi)*prec)then



         it=it+1
         ncall=ncall+inccall
         ren=dble(ncall)

         CALL VEGAS2(AVGI,SD,CHI2A)
c         CALL VEGAS2(cs,AVGI,SD,CHI2A)


         call header_out(avgi,sd)

         if(histol)then

            do j=1,nhist
               call histo2(j,0)
            enddo

         endif

         if(it.gt.itend)goto 778

         goto 777

      endif

 778  unw=.true.
      calcmax=.true.


cccccccccccccc

      if(genunw)then

         avgio=avgi
         sdo=sd

      print*,''
      print*,'**********************************************************
     &**************'
      print*,'Generating unweighted events'
      print*,'**********************************************************
     &**************'
      print*,''

c      ncall=nev
      ncall=ncallu
      if(ncall.lt.1000)ncall=1000
 566  ren=dble(ncall)
      itmx=1

      CALL VEGAS2(AVGI,SD,CHI2A)
c      CALL VEGAS2(cs,AVGI,SD,CHI2A)

      if(evnum.lt.nev)then

      print*,''
      print*,'**********************************************************
     &**************'
      write(6,100)evnum,nev
      print*,'**********************************************************
     &**************'

      endif

 100  format('  generated events so far = ',i7,'    total = ',i7)

c      ncall=ncall+nev
      ncall=ncall+inccall

      if(evnum.lt.nev)goto 566

      xsecup(1)=avgi
      xerrup(1)=sd
      xmaxup(1)=avgi

      if(enew)then
         call unwprintq
      else
         call unwprint
      endif

      if(readwt)then
      elseif((sdo/avgio).lt.(sd/avgi))then
         avgi=avgio
         sd=sdo
      endif

      endif

      close(55)
      close(45)

cccccccccccccccc

      call header_out(avgi,sd)

      if(histol)then

      do j=1,nhist
           call histo2(j,0)
      enddo

      endif

      print*,''
      write(6,301)avgi,sd
      print*,''

 301  format(' Cross section = ',G16.7,' +/-',G16.7,' pb')

      call cpu_time(t2)
      print*,'time elapsed = ', t2, ' s'
      end
