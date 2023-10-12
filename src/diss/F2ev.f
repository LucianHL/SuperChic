      subroutine F1F2ev(x,q2,f1,f2)
      implicit double precision(a-y)
      double precision garr(-6:6)

cccccc lo test
      
      call evolvePDFphoton(x,dsqrt(q2),garr,photdis1)
      
      sq=(garr(2)+garr(-2)+garr(4)+garr(-4))*4d0+garr(1)+
     &     garr(-1)+garr(3)+garr(-3)
      sq=(sq+garr(5)+garr(-5))/9d0

c      sq=(garr(2)+garr(4))*4d0/9d0
c      sq=(garr(1)+garr(3)+garr(5))/9d0
      
      f2=sq
      
c      f1=0d0
      f1=f2/2d0/x

      
      return
      end

      subroutine F1F2ap_old(x,q2,f2,fl)
      implicit double precision(a-y)

      f2=StructureFunctionxQ("EM","F2","total",x,dsqrt(q2))
      fl=StructureFunctionxQ("EM","FL","total",x,dsqrt(q2))
      
      return
      end

      subroutine F1F2ap(x,q2,f2,fl)
      implicit double precision(a-y)
      integer ie,i
      double precision garr(-6:6)

      include 'pdfinf.f'

      call lhapdf_xfxq2(2,0,900,x,q2,f2)
      call lhapdf_xfxq2(2,0,901,x,q2,fl)

      return
      end
