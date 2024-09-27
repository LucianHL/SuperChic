
ccc   prints out header information
      subroutine header_out(avgi,sd)
      implicit none
      integer outl
      double precision avgi,sd

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
      include 'beam.f'
      include 'rech.f'
      include 'diff.f'
      include 'quarkonia.f'

      call length(outtag,outl)
      open(56,file='outputs/output'//outtag(1:outl)//'.dat')
      write(56,*)''
      write(56,*)'******************************************************
     &**************'
      call length(procn,outl)
      write(56,*)'* ',procn(1:outl)
      write(56,*)'******************************************************
     &**************'
      write(56,*)
      write(56,301)avgi,sd
      write(56,*)
      write(56,*)'********************  Input parameters ***************
     &**************'
      write(56,99)' *',rts,' :  CMS collision energy (GeV)'
      write(56,97)' *',isurv,' :  Model of soft survival'
      call length(PDFname,outl)
      write(56,96)' *',PDFname(1:outl),' :  PDF set'
      write(56,97)' *',PDFmember,' :  PDF member'
      write(56,97)' *',proc,' :  Process number'
      call length(outtag,outl)
      write(56,96)' *',outtag(1:outl)//'.dat',' :  Output file'
      write(56,98)' *',sfaci,' :  Include soft survival effects'
      write(56,*)'******************************************************
     &**************'
      write(56,95)' *',diff,' : dissociation flag'
      write(56,*)'******************************************************
     &**************'
      write(56,*)''
      write(56,*)'****************** Integration parameters  ***********
     &**************'
      write(56,97)' *',ncall,' :  Preconditioning calls'
      write(56,97)' *',itmx,' :  Preconditioning iterations'
      write(56,99)' *',prec*100d0,' :  Percentage accuracy'
      write(56,97)' *',ncall1,' :  Calls in first main iteration'
      write(56,97)' *',inccall,' :  Increase calls per iteration'
      write(56,97)' *',itend,' :  Maximum number of iterations'
      write(56,97)' *',iseed,' :  Random number seed'
      write(56,*)'******************************************************
     &**************'
      write(56,97)' *',s2int,' :  Survival factor integration par.'
      write(56,*)'******************************************************
     &**************'
      write(56,*)''
      write(56,*)'********************* Unweighted Events  *************
     &**************'
      write(56,98)' *',genunw,' :  Generate unweighted events'
      if(genunw)then
         call length(erec,outl)
         write(56,97)' *',nev,' :  Number of events'
         write(56,96)' *',erec(1:outl),' :  Record format'
      endif
      write(56,*)'******************************************************
     &**************'
      write(56,*)''
      write(56,*)'********************* General Cuts *******************
     &**************'
      write(56,99)' *',ymin,' :  Minimum object rapidity'
      write(56,99)' *',ymax,' :  Maximum object rapidity'
      write(56,99)' *',mmin,' :  Minimum object mass'
      write(56,99)' *',mmax,' :  Maximum object mass'
      write(56,98)' *',gencuts,' :  Include further cuts'
      write(56,98)' *',gencuts,' :  Include spin correlations'
      write(56,*)'******************************************************
     &**************'
      write(56,*)''
      if(dps.eq.2.or.decay2)then
      write(56,*)'***************** 2-body decay cuts ******************
     &**************'
      write(56,99)' *',ptamin,' :  pT(a) min'
      write(56,99)' *',ptbmin,' :  pT(b) min'
      write(56,99)' *',etaamin,' :  eta(a) min'
      write(56,99)' *',etaamax,' :  eta(a) max'
      write(56,99)' *',etabmin,' :  eta(b) min'
      write(56,99)' *',etabmax,' :  eta(b) max'
      write(56,*)'******************************************************
     &**************'
      write(56,*)''
      elseif(dps.eq.3.or.decay3)then
      write(56,*)'***************** 3-body decay cuts ******************
     &**************'
      write(56,99)' *',ptamin3,' :  pT(a) min'
      write(56,99)' *',ptbmin3,' :  pT(b) min'
      write(56,99)' *',ptcmin3,' :  pT(c) min'
      write(56,99)' *',etaamin3,' :  eta(a) min'
      write(56,99)' *',etaamax3,' :  eta(a) max'
      write(56,99)' *',etabmin3,' :  eta(b) min'
      write(56,99)' *',etabmax3,' :  eta(b) max'
      write(56,99)' *',etacmin3,' :  eta(c) min'
      write(56,99)' *',etacmax3,' :  eta(c) max'
      write(56,*)'******************************************************
     &**************'
      write(56,*)''
      write(56,*)''
      elseif(decay4)then
      write(56,*)'***************** 4-body decay cuts ******************
     &**************'
      write(56,99)' *',ptamin4,' :  pT(a) min'
      write(56,99)' *',ptbmin4,' :  pT(b) min'
      write(56,99)' *',ptcmin4,' :  pT(c) min'
      write(56,99)' *',ptdmin4,' :  pT(d) min'
      write(56,99)' *',etaamin4,' :  eta(a) min'
      write(56,99)' *',etaamax4,' :  eta(a) max'
      write(56,99)' *',etabmin4,' :  eta(b) min'
      write(56,99)' *',etabmax4,' :  eta(b) max'
      write(56,99)' *',etacmin4,' :  eta(c) min'
      write(56,99)' *',etacmax4,' :  eta(c) max'
      write(56,99)' *',etadmin4,' :  eta(d) min'
      write(56,99)' *',etadmax4,' :  eta(d) max'
      write(56,*)'******************************************************
     &**************'
      write(56,*)''
      endif
      if(proc.eq.5.or.proc.eq.6)then
       write(56,*)'*********************** Jet cuts ********************
     &**************'
      write(56,99)' *',rjet,' :  Jet Radius'
      call length(jalg,outl)
      write(56,96)' *',jalg(1:outl),' : Record format'
      write(56,*)'******************************************************
     &**************'
      write(56,*)''
      endif
      if(proc.gt.40.and.proc.lt.46)then
      write(56,*)'**************** chi_b 2-body decays *****************
     &**************'
      write(56,99)' *',m2b,' :  Decay particle mass'
      write(56,97)' *',pdgid(6),' :  PDG number particle 1'
      write(56,97)' *',pdgid(7),' :  PDG number particle 2'
      write(56,*)'******************************************************
     &**************'
      write(56,*)''
      endif
      if(proc.gt.22.and.proc.lt.28)then
      write(56,*)'**************** chi_c 2-body decays *****************
     &**************'
      write(56,99)' *',m2b,' :  Decay particle mass'
      write(56,97)' *',pdgid(6),' :  PDG number particle 1'
      write(56,97)' *',pdgid(7),' :  PDG number particle 2'
      write(56,*)'******************************************************
     &**************'
      write(56,*)''
      endif

 99   format(a,f24.4,8x,a)
 97   format(a,i24,8x,a)
 96   format(a,a24,8x,a)
 98   format(a,l24,8x,a)
 95   format(a,a,8x,a)
 301  format(' Cross section = ',G16.7,' +/-',G16.7,' pb')

      close(56)

      return
      end
