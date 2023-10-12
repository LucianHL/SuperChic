ccc   gg --> chi_0 subprocess amplitude
      subroutine chi0(mx,mqq,qt1,qt2,out)
      implicit none
      double precision qt1sq,qt2sq,q1q2
      double precision cchi,cpp,mx,mqq
      double precision qt1(2),qt2(2)
      complex*16 out

      include 'pi.f'
      include 'quarkonia.f'

      qt1sq=qt1(1)**2+qt1(2)**2
      qt2sq=qt2(1)**2+qt2(2)**2
      q1q2=(mx**2+qt1sq+qt2sq)/2d0
  
      cchi=dsqrt(pi*mx**3*gamchi0/3d0)
      
      cpp=-(qt1(1)*qt2(1)+qt1(2)*qt2(2))*(3d0*mx**2+qt1sq+qt2sq)
      cpp=cpp-2d0*qt1sq*qt2sq
      cpp=cpp/dsqrt(6d0)*mx/2d0*cchi
      cpp=cpp/(2d0*mqq*mx+qt1sq+qt2sq)**2*4d0
      cpp=cpp*4d0*mqq**2/mx**2

      out=cpp

      return
      end
