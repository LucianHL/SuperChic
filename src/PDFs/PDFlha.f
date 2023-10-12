ccc   gluon PDF for PDFINPUT = LHAPDF
      function xgi(x,qsq)
      implicit double precision (a-z)
      double precision garr(-6:6)

      call evolvePDF(x,dsqrt(qsq),garr)

      xgi=garr(0)
 
      return
      end
