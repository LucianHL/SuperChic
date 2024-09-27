      subroutine inpdf

      include 'pdfinf.f'

cccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccc   Init LHAPDF
cccccccccccccccccccccccccccccccccccccccccccccccccccc



      call initpdfset
     &     ('PDFsets/'//PDFname)
      call initpdf(PDFmember)
      call getq2min(PDFmember,pdfq2min)
      pdfwarn=10
      return
      end
