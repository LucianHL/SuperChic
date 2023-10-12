ccc   strong coupling constant, for PDFINPUT = LHAPDF
      function alphas(qsq)
      implicit double precision (a-z)

      alphas=alphasPDF(dsqrt(qsq))

      return
      end
