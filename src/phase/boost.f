ccc   boosting subroutine
      subroutine boost(p,pboo,pcm,plb)
      real*8  pboo(4),pcm(4),plb(4),p,fact
         plb(4)=(pboo(4)*pcm(4)+pboo(3)*pcm(3)
     &             +pboo(2)*pcm(2)+pboo(1)*pcm(1))/p
         fact=(plb(4)+pcm(4))/(p+pboo(4))
         do 10 j=1,3
  10     plb(j)=pcm(j)+fact*pboo(j)
      return
      end
