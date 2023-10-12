ccc   gg --> eta_{c,b} subprocess amplitude
      subroutine etaq(mx,qt1,qt2,out)
      implicit none
      double precision qt1sq,qt2sq,q1q2,cpp,ceta,mx
      double precision qt1(2),qt2(2)
      complex*16 out

      include 'pi.f'
      include 'zi.f'
      include 'quarkonia.f'

      qt1sq=qt1(1)**2+qt1(2)**2
      qt2sq=qt2(1)**2+qt2(2)**2
      q1q2=(mx**2+qt1sq+qt2sq)/2d0

      ceta=dsqrt(2d0*pi*mx*gameta)
      
      cpp=(qt1(1)*qt2(2)-qt1(2)*qt2(1))
      cpp=cpp*mx**2/2d0*ceta
      cpp=cpp/q1q2

ccccccccccccc

      out=cpp*zi

      return
      end
