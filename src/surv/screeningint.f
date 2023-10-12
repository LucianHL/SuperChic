ccc   interpolator for screened amplitude
      subroutine screeningint(i,j,ktsq,out,out1)
      implicit none
      double precision m,m1,ktmin,inckt,del1,del
      double precision out,out1,ktsq
      integer i,j,it

      include 'screenamp.f'

      if(ktsq.lt.0.001d0)then

         it=0

         m=(sca(i,j,it+1,2)-sca(i,j,it,2))
     &        /(dexp(sca(i,j,it+1,1))-sca(i,j,it,1))
         del=ktsq-sca(1,1,it,1)
         m1=sca1(i,j,it+1,2)-sca1(i,j,it,2)
     &        /(dexp(sca1(i,j,it+1,1))-sca1(i,j,it,1))
         del1=ktsq-sca1(1,1,it,1)

         out=m*del+sca(i,j,it,2)
         out1=m1*del1+sca1(i,j,it,2)

      elseif(ktsq.lt.8d0)then

         ktmin=sca(1,1,1,1)
         inckt=sca(1,1,2,1)-sca(1,1,1,1)
         it=nint((dlog(ktsq)-ktmin)/inckt)
         if(dble(it).gt.((dlog(ktsq)-ktmin)/inckt))then
            it=it-1
         endif

         m=(sca(i,j,it+2,2)-sca(i,j,it+1,2))
     &        /(sca(i,j,it+2,1)-sca(i,j,it+1,1))
         del=dlog(ktsq)-sca(1,1,it+1,1)
         m1=sca1(i,j,it+2,2)-sca1(i,j,it+1,2)
     &        /(sca1(i,j,it+2,1)-sca1(i,j,it+1,1))
         del1=dlog(ktsq)-sca1(1,1,it+1,1)

         out=m*del+sca(i,j,it+1,2)
         out1=m1*del1+sca1(i,j,it+1,2)

      else

         out=0d0
         out1=0d0

      endif

      return
      end
