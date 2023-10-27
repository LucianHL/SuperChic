ccc   strong coupling constant, for PDFINPUT = LHAPDF
      function alphas(qsq)
      implicit double precision (a-z)
      if (dsq.le.100.d0) then 
      alphas=alphasPDF(10.0d0)
      else
      alphas=alphasPDF(dsqrt(qsq))
      end if
      return
      end
