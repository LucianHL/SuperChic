      function betaion(t)
      implicit none
      double precision t,a,b,c,expon,betaion

      a=0.167d0
      b=22d0
      c=0.402d0

      a=0.18d0
      b=20.8d0
      c=0.414d0

      expon=(b*(a-t))**c-(a*b)**c
      betaion=dexp(-expon)

      return
      end
