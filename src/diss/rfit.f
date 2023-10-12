      function rfit(i,x,qsq)
      implicit double precision(a-z)
      integer i

      theta=1d0+12d0*(qsq/(qsq+1d0))*(0.125d0**2/(0.125d0**2+x**2))
      
      if(i.eq.1)then
         a1=0.0485d0
         a2=0.5470d0
         a3=2.0621d0
         a4=-0.3804d0
         a5=0.5090d0
         a6=-0.0285d0

         rfit=a1/dlog(qsq/0.04d0)*theta
         rfit=rfit+a2*(1d0+a4*x+a5*x**2)*x**a6/(qsq**4+a3**4)**0.25d0
         
      elseif(i.eq.2)then
         b1=0.0481d0
         b2=0.6114d0
         b3=-0.3509d0
         b4=-0.4611d0
         b5=0.7172d0
         b6=-0.0317d0

         rfit=b1/dlog(qsq/0.04d0)*theta
         rfit=rfit+(b2/qsq+b3/(qsq**2+0.3d0**2))*
     &        (1d0+b4*x+b5*x**2)*x**b6
         
      elseif(i.eq.3)then
         c1=0.0577d0
         c2=0.4644d0
         c3=1.8288d0
         c4=12.3708d0
         c5=-43.1043d0
         c6=41.7415d0

         qthr=c4*x+c5*x**2+c6*x**3
         
         rfit=c1/dlog(qsq/0.04d0)*theta
         rfit=rfit+c2/dsqrt((qsq-qthr)**2+c3**2)

      endif

      return
      end
