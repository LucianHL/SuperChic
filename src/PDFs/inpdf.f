      subroutine inpdf

      include 'pdfinf.f'

cccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccc   Init LHAPDF
cccccccccccccccccccccccccccccccccccccccccccccccccccc

  

      call initpdfset
     &     ('PDFsets/'//PDFname)
      call initpdf(PDFmember)

      return
      end
