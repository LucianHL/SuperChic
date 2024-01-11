ccc   strong coupling constant, for PDFINPUT = LHAPDF
      function alphas(qsq)
      implicit double precision (a-z)
      include 'pdfinf.f'

      if (qsq.lt.pdfq2min) then 
      if ( pdfwarn.ge.0) then
ccc   LHL 11/1/24 omit these to avoid confusion for now
c      write(*,*)'Warning in alphas.f qsq = ',qsq,' < ',pdfq2min,pdfwarn
c      write(*,*)'alphas being frozen at qsq = ',pdfq2min
      pdfwarn=pdfwarn-1
      endif
      alphas=alphasPDF(pdfq2min)
      else
      alphas=alphasPDF(dsqrt(qsq))
      end if

      return
      end
