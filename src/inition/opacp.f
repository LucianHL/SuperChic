      function opacp(btt)
      implicit double precision(a-y)
      integer n,ntot

      include 'pi.f'
      include 'onechannel.f'
     
      qtmax=2d0
      ntot=1000
      hq=qtmax/dble(ntot)

      sum=0d0
      
      do n=1,ntot

         qt=(dble(n)-0.5d0)*hq

         wt=qt*hq/2d0/pi
         wt=wt*betaion(-qt**2)**2*besj0(qt*btt)*sig0
         
         sum=sum+wt
         
      enddo

      opacp=sum
      
      return
      end
