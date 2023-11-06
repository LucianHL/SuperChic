ccc   prints out header information
      subroutine headerlhe
      implicit double precision(a-y)
      integer outl

      include 'gencuts.f'
      include 'range.f'
      include 'unweighted.f'
      include 'vegas.f'
      include 'vegaspars.f'
      include 'proc.f'
      include 'pdfinf.f'
      include 'genunw.f'
      include 'decay.f'
      include 'survin.f'
      include 'procn.f'
      include 'vars.f'
      include 'pdg.f'
      include 'jetalg.f'
      include 'record.f'
      include 'output.f'
      include 'nsurv.f'
      include 'prec.f'
      include 'intag.f'

      call length(procn,outl)

      write(45,*)'***********************************************'
      write(45,*)'************** Superchic v4.2 *****************'
      write(45,*)'***********************************************'
      write(45,*)'*  v4.2               DATE    13/03/23        *'
      write(45,*)'*                                             *'
      write(45,*)'*  Author: Lucian Harland-Lang                *'
      write(45,*)'*  (l.harland-lang@ucl.ac.uk)                 *'
      write(45,*)'*                                             *'
      write(45,*)'*  For details see :                          *'
      write(45,*)'*                                             *'
      write(45,*)'*  arXiv 2303.04826 (ion dissociation)        *'
      write(45,*)'*  arXiv 2201.08403 (WW updates)              *'
      write(45,*)'*  arXiv 2007.12704 (v4 updates)              *'
      write(45,*)'*  arXiv 1812.04886 (SUSY)                    *'
      write(45,*)'*  arXiv 1810.06567                           *'
      write(45,*)'*  arXiv 1508.02718                           *'
      write(45,*)'*  arXiv 1405.0018 (review)                   *'
      write(45,*)'*  arXiv 1005.0695 (quarkonia and diphoton    *'
      write(45,*)'*  arXiv 1011.0680 (quarkonia)                *'
      write(45,*)'*  arXiv 1105.1626 (meson pairs)              *'
      write(45,*)"*  arXiv 1302.2004 (meson pairs - eta/eta')   *"
      write(45,*)'*  arXiv 1306.6661 (Skewed PDF)               *'
      write(45,*)'*  arXiv 1306.2149 (Skewed PDF)               *'
      write(45,*)'*                                             *'
      write(45,*)'*  Available at :                             *'
      write(45,*)'*  http:://projects.hepforge.org/superchic    *'
      write(45,*)'*                                             *'
      write(45,*)'***********************************************'
      write(45,*)''
      write(45,*)'******************************************************
     &******************'
      call length(procn,outl)
      write(45,*)'* ',procn(1:outl),' production'
      write(45,*)'******************************************************
     &******************'
      write(45,*)''
      write(45,*)'********************  Input parameters ***************
     &******************'
      write(45,99)' *',rts,' :  CMS collision energy (GeV)'
      write(45,97)' *',isurv,' :  Model of soft survival'
      call length(intag,outl)
      write(45,96)' *',intag(1:outl),' :  input file'
      call length(PDFname,outl)
      write(45,96)' *',PDFname(1:outl),' :  PDF set'
      write(45,97)' *',PDFmember,' :  PDF member'
      write(45,97)' *',proc,' :  Process number'
      call length(outtag,outl)
      write(45,96)' *',outtag(1:outl)//'.dat',' :  Output file'
      write(45,98)' *',sfaci,' :  Include soft survival effects'
      write(45,*)'******************************************************
     &******************'
      write(45,*)''
      write(45,*)'********************* Integration parameters  ********
     &******************'
      write(45,97)' *',ncall,' :  Preconditioning calls'
      write(45,97)' *',itmx,' :  Preconditioning iterations'
      write(45,99)' *',prec,' :  Percentage accuracy'
      write(45,97)' *',ncall1,' :  Calls in first main iteration'
      write(45,97)' *',inccall,' :  Increase calls per iteration'
      write(45,97)' *',itend,' :  Maximum number of iterations'
      write(45,97)' *',iseed,' :  Random number seed'
      write(45,*)'******************************************************
     &******************'
      write(45,97)' *',s2int,' : Survival factor integration par.'
      write(45,*)'******************************************************
     &******************'
      write(45,*)''
      write(45,*)'********************* Unweighted Events  *************
     &******************'
      write(45,98)' *',genunw,' :  Generate unweighted events'
      if(genunw)then
         call length(erec,outl)
         write(45,97)' *',nev,' :  Number of events'
         write(45,96)' *',erec(1:outl),' : Record format'
      endif
      write(45,*)'******************************************************
     &******************'
      write(45,*)''
      write(45,*)'********************* General Cuts *******************
     &******************'
      write(45,99)' *',ymin,' :  Minimum object rapidity'
      write(45,99)' *',ymax,' :  Maximum object rapidity'
      write(45,99)' *',mmin,' :  Minimum object mass'
      write(45,99)' *',mmax,' :  Maximum object mass'
      write(45,98)' *',gencuts,' :  Include further cuts'
      write(45,98)' *',gencuts,' :  Include spin correlations'
      write(45,*)'*****************************************************
     &******************'
      write(45,*)''
      if(decay4)then
      write(45,*)'***************** 4-body decay cuts ******************
     &******************'
      write(45,99)' *',ptamin4,' :  pT(a) min'
      write(45,99)' *',ptbmin4,' :  pT(b) min'
      write(45,99)' *',ptcmin4,' :  pT(c) min'
      write(45,99)' *',ptdmin4,' :  pT(d) min'
      write(45,99)' *',etaamin4,' :  eta(a) min'
      write(45,99)' *',etaamax4,' :  eta(a) max'
      write(45,99)' *',etabmin4,' :  eta(b) min'
      write(45,99)' *',etabmax4,' :  eta(b) max'
      write(45,99)' *',etacmin4,' :  eta(c) min'
      write(45,99)' *',etacmax4,' :  eta(c) max'
      write(45,99)' *',etadmin4,' :  eta(d) min'
      write(45,99)' *',etadmax4,' :  eta(d) max'
      write(45,*)'******************************************************
     &******************'
      write(45,*)''
      elseif(dps.eq.2.or.decay2)then
      write(45,*)'***************** 2-body decay cuts ******************
     &*****************'
      write(45,99)' *',ptamin,' :  pT(a) min'
      write(45,99)' *',ptbmin,' :  pT(b) min'
      write(45,99)' *',etaamin,' :  eta(a) min'
      write(45,99)' *',etaamax,' :  eta(a) max'
      write(45,99)' *',etabmin,' :  eta(b) min'
      write(45,99)' *',etabmax,' :  eta(b) max'
      write(45,*)'******************************************************
     &*****************'
      write(45,*)''
      elseif(dps.eq.3.or.decay3)then
      write(45,*)'***************** 3-body decay cuts ******************
     &******************'
      write(45,99)' *',ptamin3,' :  pT(6) min'
      write(45,99)' *',ptbmin3,' :  pT(7) min'
      write(45,99)' *',ptcmin3,' :  pT(8) min'
      write(45,99)' *',etaamin3,' :  eta(a) min'
      write(45,99)' *',etaamax3,' :  eta(a) max'
      write(45,99)' *',etabmin3,' :  eta(b) min'
      write(45,99)' *',etabmax3,' :  eta(b) max'
      write(45,99)' *',etacmin3,' :  eta(c) min'
      write(45,99)' *',etacmax3,' :  eta(c) max'
      write(45,*)'******************************************************
     &******************'
      write(45,*)''
      endif
      if(proc.eq.7.or.proc.eq.8)then
       write(45,*)'*********************** Jet cuts ********************
     &******************'
      write(45,99)' *',rjet,' :  Jet Radius'
      call length(jalg,outl)
      write(45,96)' *',jalg(1:outl),' :  Jet alg.'
      write(45,*)'******************************************************
     &******************'
      write(45,*)''
      endif
      if(proc.gt.41.and.proc.lt.47)then
      write(45,*)'**************** chi_b 2-body decays *****************
     &****************'
      write(45,99)' *',m2b,' :  Decay particle mass'
      write(45,97)' *',pdgid(6),' :  PDG number particle 1'
      write(45,97)' *',pdgid(7),' :  PDG number particle 2'
      write(45,*)'******************************************************
     &******************'
      write(45,*)''
      endif
      if(proc.gt.23.and.proc.lt.29)then
      write(45,*)'**************** chi_c 2-body decays *****************
     &****************'
      write(45,99)' *',m2b,' :  Decay particle mass'
      write(45,97)' *',pdgid(6),' :  PDG number particle 1'
      write(45,97)' *',pdgid(7),' :  PDG number particle 2'
      write(45,*)'******************************************************
     &******************'
      endif

 99   format(a,f24.4,8x,a)
 97   format(a,i24,8x,a)
 96   format(a,a24,8x,a)
 98   format(a,l24,8x,a)

      return
      end
