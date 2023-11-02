ccc   gg --> chi_0 subprocess amplitude
      subroutine chi2(p,mqq,mxx,q1,q2,echi2,out)
      implicit none
      double precision qt1sq,qt2sq,q1q2,cchi,mqq,mxx
      double precision q1(2),q2(2)
      double precision qt1(4),qt2(4)
      complex*16 out,echi2(5,4,4),cpp
      integer p,i,j
      double precision pcm(4),pboo(4),plb(4)
      double precision q1b(4),q2b(4)

      include 'pi.f'
      include 'zi.f'
      include 'vars.f'
      include 'mom.f'
      include 'quarkonia.f'

      qt1(4)=0d0
      qt2(4)=0d0
      qt1(3)=0d0
      qt2(3)=0d0

      do i=1,2
         qt1(i)=q1(i)
         qt2(i)=q2(i)
      enddo

      do i=1,4
         pboo(i)=-q(i,5)
      enddo
         pboo(4)=q(4,5)

      do i=1,4
         pcm(i)=q(i,1)
      enddo
      call boost(mx,pboo,pcm,plb)
      do i=1,4
         q1b(i)=plb(i)
      enddo

      do i=1,4
         pcm(i)=q(i,2)
      enddo
      call boost(mx,pboo,pcm,plb)
      do i=1,4
         q2b(i)=plb(i)
      enddo

      call boost(mx,pboo,qt1,plb)
      do i=1,4
         qt1(i)=plb(i)
      enddo

      call boost(mx,pboo,qt2,plb)
      do i=1,4
         qt2(i)=plb(i)
      enddo

      qt1sq=-(qt1(4)**2-qt1(3)**2-qt1(2)**2-qt1(1)**2)
      qt2sq=-(qt2(4)**2-qt2(3)**2-qt2(2)**2-qt2(1)**2)
      q1q2=(mx**2+qt1sq+qt2sq)/2d0

      cchi=dsqrt(pi*mx**3*gamchi0/3d0)

      do i=1,3
         qt1(i)=-qt1(i)
         qt2(i)=-qt2(i)
         q1b(i)=-q1b(i)
         q2b(i)=-q2b(i)
      enddo

      cpp=(0d0,0d0)

      do i=1,4
         do j=1,4
            cpp=cpp+s*qt1(i)*qt2(j)*echi2(p,i,j)-2d0*(qt1(1)*qt2(1)+
     &           qt1(2)*qt2(2)+qt1(3)*qt2(3)-qt1(4)*qt2(4))*
     &           q1b(i)*q2b(j)*echi2(p,i,j)
         enddo
      enddo

      cpp=cpp*cchi*dsqrt(2d0)*mx/s
      cpp=cpp/(2d0*mqq*mx+qt1sq+qt2sq)**2*4d0
      cpp=cpp*dsqrt(mx/mxx)
      cpp=cpp*mx**2/2d0

      out=cpp

      return
      end
