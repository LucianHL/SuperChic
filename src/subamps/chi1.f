ccc   gg --> chi_1 subprocess amplitude
      subroutine chi1(p,mx,mqq,mxx,qt1,qt2,e1,out)
      implicit none
      double precision q1q2,qt1sq,qt2sq,cchi,mxx,mqq,mx
      double precision qt1(2),qt2(2)
      complex*16 out,e1(3,4),cpp
      integer p

      include 'pi.f'
      include 'zi.f'
      include 'quarkonia.f'     

      qt1sq=qt1(1)**2+qt1(2)**2
      qt2sq=qt2(1)**2+qt2(2)**2

      q1q2=(mx**2+qt1sq+qt2sq)/2d0
      cchi=dsqrt(pi*mx**3*gamchi0/3d0)
      
      cpp=e1(p,1)*(qt2(2)*qt1sq-qt1(2)*qt2sq)-
     &e1(p,2)*(qt2(1)*qt1sq-qt1(1)*qt2sq)

      cpp=cpp*cchi*zi
      cpp=cpp/(2d0*mqq*mx+qt1sq+qt2sq)**2*4d0
      cpp=cpp*4d0*mqq**2/mx**2
      cpp=cpp*dsqrt(mx/mxx)
      cpp=cpp*mx**2/2d0

      out=cpp

      return
      end
