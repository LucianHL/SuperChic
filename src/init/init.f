ccc   Initialises grids for skewed PDFs and survival factors
      implicit double precision(a-y)
      integer isurv
      character*100 dum

      include 'pi.f'
      include 'vars.f'
      include 'intag.f'
      include 'pdfinf.f'
      include 'mp.f'

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
      
c      stop
      
      call calcsud           ! sudakov factor
      call calchg            ! skewed PDF

      print*,'Now run ./superchic'

      stop
      end
