      function rhozxy(rxy)
      implicit none
      double precision sum,wt,hz,rrz,r,rhoz
      double precision rhozxy,rzmin,rzmax,rxy
      integer i,itot

      include 'ion.f'
      
      itot=10000

      rzmax=10d0*rzg
      hz=rzmax/dble(itot)
      
      sum=0d0
      
      do i=1,itot

         rrz=(dble(i)-0.5d0)*hz

         r=dsqrt(rrz**2+rxy**2)
         wt=rhoz(r)*hz
         
         sum=sum+wt

      enddo

      rhozxy=sum*2d0

      return
      end

      function rhonxy(rxy)
      implicit none
      double precision sum,wt,hn,rrn,r,rhon
      integer i,itot
      double precision rhonxy,rnmin,rnmax,rxy

      include 'ion.f'
      
      itot=10000

      rnmax=10d0*rng
      hn=rnmax/dble(itot)
      
      sum=0d0
      
      do i=1,itot

         rrn=(dble(i)-0.5d0)*hn

         r=dsqrt(rrn**2+rxy**2)
         wt=rhon(r)*hn
         
         sum=sum+wt

      enddo

      rhonxy=sum*2d0

      return
      end


      function rhozxyc(rxy)
      implicit none
      double precision sum,wt,hz,rrz,r,rhoz
      integer i,itot
      double precision rhozxyc,rzmin,rzmax,rxy

      include 'ion.f'
      
      itot=10000

      rzmax=10d0*rzg
      rzmin=rcore
      hz=(rzmax-rzmin)/dble(itot)
      
      sum=0d0
      
      do i=1,itot

         rrz=(dble(i)-0.5d0)*hz+rzmin

         r=dsqrt(rrz**2+rxy**2)
         wt=rhoz(r)*hz
         
         sum=sum+wt

      enddo

      rhozxyc=sum*2d0

      return
      end

      function rhonxyc(rxy)
      implicit none
      double precision sum,wt,hn,rrn,r,rhon
      integer i,itot
      double precision rhonxyc,rnmin,rnmax,rxy

      include 'ion.f'
      
      itot=10000

      rnmax=10d0*rng
      rnmin=rcore
      hn=(rnmax-rnmin)/dble(itot)
      
      sum=0d0
      
      do i=1,itot

         rrn=(dble(i)-0.5d0)*hn+rnmin

         r=dsqrt(rrn**2+rxy**2)
         wt=rhon(r)*hn
         
         sum=sum+wt

      enddo

      rhonxyc=sum*2d0

      return
      end
