
ccc   prints out header information
      subroutine header
      implicit none
      integer outl,outl1, i

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
      include 'head.f'

      call length(procn,outl)
      do i=1,30
      print*,head(i)
      end do
      
      call length(procn,outl)
      call length(beam,outl1)
      print*,'* ',procn(1:outl),' production'
     &     ,' (',beam(1:outl1),' beam)'
      print*,'**********************************************************
     &**************'
      print*,''
      print*,'********************  Input parameters *******************
     &**************'
      write(*,99)' *',rts,' :  CMS collision energy (GeV)'
      write(*,97)' *',isurv,' :  Model of soft survival'
      call length(intag,outl)
      write(*,96)' *',intag(1:outl),' :  input file'
      call length(PDFname,outl)
      write(*,96)' *',PDFname(1:outl),' :  PDF set'
      write(*,97)' *',PDFmember,' :  PDF member'
      write(*,97)' *',proc,' :  Process number'
      call length(outtag,outl)
      write(*,96)' *',outtag(1:outl)//'.dat',' :  Output file'
      write(*,98)' *',sfaci,' :  Include soft survival effects'
      print*,'**********************************************************
     &**************'
      if(diffsd.eq.'n')then
         write(*,96)' *',diff(1:3),' : Dissociation flag'
      else
         write(*,96)' *',diffsd(1:3),' : Dissociation flag'
      endif
      print*,'**********************************************************
     &**************'
      print*,''
      print*,'********************* Integration parameters  ************
     &**************'
      write(*,97)' *',ncall,' :  Preconditioning calls'
      write(*,97)' *',itmx,' :  Preconditioning iterations'
      write(*,99)' *',prec,' :  Percentage accuracy'
      write(*,97)' *',ncall1,' :  Calls in first main iteration'
      write(*,97)' *',inccall,' :  Increase calls per iteration'
      write(*,97)' *',itend,' :  Maximum number of iterations'
      write(*,97)' *',iseed,' :  Random number seed'
      print*,'**********************************************************
     &**************'
      write(*,97)' *',s2int,' : Survival factor integration par.'
      print*,'**********************************************************
     &**************'
      print*,''
      print*,'********************* Unweighted Events  *****************
     &**************'
      write(*,98)' *',genunw,' :  Generate unweighted events'
      if(genunw)then
         call length(erec,outl)
         write(*,97)' *',nev,' :  Number of events'
         if(erech)then
            write(*,96)' *','hepmc',' : Record format'
         else
            write(*,96)' *',erec(1:outl),' : Record format'
         endif
      endif
      print*,'**********************************************************
     &**************'
      print*,''
      print*,'********************* General Cuts ***********************
     &**************'
      write(*,99)' *',ymin,' :  Minimum object rapidity'
      write(*,99)' *',ymax,' :  Maximum object rapidity'
      write(*,99)' *',mmin,' :  Minimum object mass'
      write(*,99)' *',mmax,' :  Maximum object mass'
      write(*,98)' *',gencuts,' :  Include further cuts'
      write(*,98)' *',gencuts,' :  Include spin correlations'
      print*,'**********************************************************
     &**************'
      print*,''
      if(decay4)then
      print*,'***************** 4-body decay cuts **********************
     &**************'
      write(*,99)' *',ptamin4,' :  pT(a) min'
      write(*,99)' *',ptbmin4,' :  pT(b) min'
      write(*,99)' *',ptcmin4,' :  pT(c) min'
      write(*,99)' *',ptdmin4,' :  pT(d) min'
      write(*,99)' *',etaamin4,' :  eta(a) min'
      write(*,99)' *',etaamax4,' :  eta(a) max'
      write(*,99)' *',etabmin4,' :  eta(b) min'
      write(*,99)' *',etabmax4,' :  eta(b) max'
      write(*,99)' *',etacmin4,' :  eta(c) min'
      write(*,99)' *',etacmax4,' :  eta(c) max'
      write(*,99)' *',etadmin4,' :  eta(d) min'
      write(*,99)' *',etadmax4,' :  eta(d) max'
      print*,'**********************************************************
     &**************'
      print*,''
      elseif(dps.eq.2.or.decay2)then
      print*,'***************** 2-body decay cuts **********************
     &**************'
      write(*,99)' *',ptamin,' :  pT(a) min'
      write(*,99)' *',ptbmin,' :  pT(b) min'
      write(*,99)' *',etaamin,' :  eta(a) min'
      write(*,99)' *',etaamax,' :  eta(a) max'
      write(*,99)' *',etabmin,' :  eta(b) min'
      write(*,99)' *',etabmax,' :  eta(b) max'
      print*,'**********************************************************
     &**************'
      print*,''
      elseif(dps.eq.3.or.decay3)then
      print*,'***************** 3-body decay cuts **********************
     &**************'
      write(*,99)' *',ptamin3,' :  pT(6) min'
      write(*,99)' *',ptbmin3,' :  pT(7) min'
      write(*,99)' *',ptcmin3,' :  pT(8) min'
      write(*,99)' *',etaamin3,' :  eta(a) min'
      write(*,99)' *',etaamax3,' :  eta(a) max'
      write(*,99)' *',etabmin3,' :  eta(b) min'
      write(*,99)' *',etabmax3,' :  eta(b) max'
      write(*,99)' *',etacmin3,' :  eta(c) min'
      write(*,99)' *',etacmax3,' :  eta(c) max'
      print*,'**********************************************************
     &**************'
      print*,''
      endif
      if(proc.eq.7.or.proc.eq.8)then
       print*,'*********************** Jet cuts ************************
     &**************'
      write(*,99)' *',rjet,' :  Jet Radius'
      call length(jalg,outl)
      write(*,96)' *',jalg(1:outl),' :  Jet alg.'
      print*,'**********************************************************
     &**************'
      print*,''
      endif
      if(proc.gt.41.and.proc.lt.47)then
      print*,'**************** chi_b 2-body decays *******************
     &**************'
      write(*,99)' *',m2b,' :  Decay particle mass'
      write(*,97)' *',pdgid(6),' :  PDG number particle 1'
      write(*,97)' *',pdgid(7),' :  PDG number particle 2'
      print*,'**********************************************************
     &**************'
      print*,''
      endif
      if(proc.gt.23.and.proc.lt.29)then
      print*,'**************** chi_c 2-body decays *******************
     &**************'
      write(*,99)' *',m2b,' :  Decay particle mass'
      write(*,97)' *',pdgid(6),' :  PDG number particle 1'
      write(*,97)' *',pdgid(7),' :  PDG number particle 2'
      print*,'**********************************************************
     &**************'
      print*,''
      endif

 99   format(a,f24.4,8x,a)
 97   format(a,i24,8x,a)
 96   format(a,a24,8x,a)
 98   format(a,l24,8x,a)

      return
      end
