c     alpha_EM(Q^2)

      function alphaEM(qsq)
      implicit none
      double precision alphaEM,qsq
      double precision alem,alem3pi,rpigg

      alem=1./137.
      alem3pi=alem/(3.0*3.1415926535)      
      if (qsq.lt.2.0e-6) then
         rpigg=0.0
      elseif (qsq.lt.0.09) then
         rpigg = alem3pi*(13.4916 + dlog(qsq)) + 0.00835*dlog(1.0 + qsq)
      elseif ( qsq.lt.9.0 ) then
         rpigg=alem3pi*(16.3200 + 2.0*dlog(qsq)) + 0.00238*dlog(1.0 +
     $        3.927*qsq)
      elseif ( qsq.lt.10000.0 ) then
         rpigg=alem3pi*(13.4955 + 3.0*dlog(qsq)) + 0.00165 + 0.00299
     $        *dlog(1.0 + qsq)
      else
         rpigg=alem3pi*(13.4955 + 3.0*dlog(qsq)) + 0.00221 + 0.00293
     $        *dlog(1.0 + qsq)
      endif
      alphaEM=alem/(1.0-rpigg)


c      alphaEM=1d0/137d0

      
      return
      end
