ccc   generates two-body W decay, with additional information
ccc   saved for spin correlations
      subroutine twobodyw(in,i1,i2,m1,m2)
      implicit none
      double precision m1,m2,ran2,ran9,ran10
      double precision sphi,sint,phi,pcmod,mx,cost,cphi
      integer in,i1,i2
      integer k
      double precision pcm(4),pboo(4),plb(4)

      include 'mom.f'
      include 'partonmom4.f'
      include 'pi.f'
      include 'elcollw.f'
 
      mx=dsqrt(q(4,in)**2-q(3,in)**2-q(2,in)**2-q(1,in)**2)

      if(elcollw)then
         if(in.eq.6)then
            ran9=r9a
            ran10=r10a
         elseif(in.eq.7)then
            ran9=r9b
            ran10=r10b
         endif
      else
         ran9=ran2()
         ran10=ran2()
         if(in.eq.6)then
            r9a=ran9
            r10a=ran10
         elseif(in.eq.7)then
            r9b=ran9
            r10b=ran10
         endif
      endif

c      print*,in,ran9,ran10
      
      cost=2d0*ran9-1d0
      phi=2d0*pi*ran10
      sint=dsqrt(1d0-cost**2)
      cphi=dcos(phi)
      sphi=dsin(phi)
      pcm(4)=(mx**2+m1**2-m2**2)/(2d0*mx)
      pcmod=dsqrt(pcm(4)**2-m1**2)
      pcm(1)=pcmod*sint*sphi
      pcm(2)=pcmod*sint*cphi
      pcm(3)=pcmod*cost   

      do k=1,4
      pboo(k)=q(k,in)
      enddo
      call boost(mx,pboo,pcm,plb)
      do k=1,4
      q(k,i1)=plb(k)
      q(k,i2)=q(k,in)-q(k,i1)
      enddo

      if(in.eq.6)then
         do k=1,4
            paa(k)=pcm(k)
         enddo
      elseif(in.eq.7)then
         do k=1,4
            pbb(k)=pcm(k)
         enddo
      endif

      return
      end
