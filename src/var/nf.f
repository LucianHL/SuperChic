ccc   sets number of active flavours
      function nf(qsq)
      implicit none
      double precision qsq,mc,mb,nf

      mc=1.4d0
      mb=4.75d0
      
      if(qsq.lt.mc**2)then
      nf=3d0
      elseif(qsq.lt.mb**2)then
      nf=4d0
      else
      nf=5d0
      endif

      return
      end
