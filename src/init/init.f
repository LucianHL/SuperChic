ccc   Initialises grids for skewed PDFs and survival factors
      program initpdfs
      implicit none
      integer isurv
      character*100 dum

      include 'pi.f'
      include 'vars.f'
      include 'intag.f'
      include 'pdfinf.f'
      include 'mp.f'
      include 'beam.f'
      include 'proc.f'

      call system('mkdir -p inputs evrecs outputs')
      mp=0.938272046d0
      pi=dacos(-1d0)

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
      read(*,*)dum
      read(*,*)proc
      read(*,*)beam
      if (beam .eq. 'el') then 
      write(*,*)'Running the init program is not required for ee beams'
      goto 999
      end if

cccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccc   Init LHAPDF
cccccccccccccccccccccccccccccccccccccccccccccccccccc

      call inpdf

cccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccc

      call initpars(isurv)   ! Initialise soft survival parameters
      call calcop               ! proton opacity
      call calcscreen        ! screening amplitude

      call readcoh
      call dd
      call readscreen
      call sdcoh
      call apfelinit
      call sdincoh


      call calcsud           ! sudakov factor
      call calchg            ! skewed PDF

 999  continue
      print*,'Now run ./superchic'

      end
