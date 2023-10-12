ccc   outputs exclusive skewed pdf
      function fg(x,qsq,mx)
      implicit double precision (a-z)

      call sudint(qsq,mx,tg,dtg)
      call hpdfint(dlog(x),qsq,hgo,diffhg)
      
      fg=dsqrt(tg)*diffhg+hgo*dtg/2d0/dsqrt(tg)
      fg=fg*qsq

  
      return
      end

     
