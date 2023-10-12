      function betaionex(t)
      implicit none
      double precision t,expon,betaionex

      include 'onechannel.f'

      expon=(b*(a-t))**c-(a*b)**c
      
      betaionex=dexp(-expon)
      
      return
      end
