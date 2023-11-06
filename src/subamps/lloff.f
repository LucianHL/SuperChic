ccc   gamma gamma --> l+l- subprocess amplitude - off-shell
      subroutine lloff(p)
      implicit none
      double precision q1(4),q2(4),pel(4),ppos(4)
      double precision q1m(4),q2m(4),pelm(4),pposm(4)
      double precision ul,tl,qsq1,qsq2,alphaem,beta,al
      complex*16 uelb(4),vpos(4),zt2,zt1,zout1,zout2,zout
      integer i1,i2
      integer i,j,l,p

      include 'mom.f'
      include 'gmatrices.f'
      include 'gmatrices_comb.f'
      include 'ewpars.f'
      include 'pi.f'
      include 'vars.f'
      include 'mandelstam.f'
      include 'partonmom2.f'
      include 'norm.f'
      include 'mq.f'
      include 'zi.f'
      include 'zoutarr.f'
      include 'eff.f'

      do i=1,4
         q1(i)=q(i,1)-q(i,3)
         q2(i)=q(i,2)-q(i,4)
         pel(i)=q(i,6)
         ppos(i)=q(i,7)
         q1m(i)=-q1(i)
         q2m(i)=-q2(i)
         pelm(i)=-q(i,6)
         pposm(i)=-q(i,7)
      enddo
      pelm(4)=q(4,6)
      pposm(4)=q(4,7)
      q1m(4)=q1(4)
      q2m(4)=q2(4)

      al=0d0

      beta=dsqrt(1d0-4d0*mq**2/mx**2)

      qsq1=(q(4,3)-q(4,1))**2-(q(3,3)-q(3,1))**2-(q(2,3)-q(2,1))**2
     &     -(q(1,3)-q(1,1))**2
      qsq1=-qsq1
      qsq2=(q(4,4)-q(4,2))**2-(q(3,4)-q(3,2))**2-(q(2,4)-q(2,2))**2
     &     -(q(1,4)-q(1,2))**2
      qsq2=-qsq2

      if(p.eq.1)then ! ++
         call upb(6,uelb)
         call vp(7,vpos)
      elseif(p.eq.2)then ! --
         call umb(6,uelb)
         call vm(7,vpos)
      elseif(p.eq.3)then ! +-
         call umb(6,uelb)
         call vp(7,vpos)
      elseif(p.eq.4)then ! -+
         call upb(6,uelb)
         call vm(7,vpos)
      endif

      ul=(q1(4)-ppos(4))**2-(q1(3)-ppos(3))**2
     &     -(q1(2)-ppos(2))**2-(q1(1)-ppos(1))**2
      tl=(q2(4)-ppos(4))**2-(q2(3)-ppos(3))**2
     &     -(q2(2)-ppos(2))**2-(q2(1)-ppos(1))**2

      do i1=1,4
         do i2=1,4

            zout1=0d0
            zout2=0d0

      do 900 i=1,4
      do 900 l=1,4

         do j=1,4

            zt1=uelb(i)*gmatrix_3(i2,j,i1,i,l)*
     &           vpos(l)*(q1(j)-ppos(j))
            if(j.lt.4)zt1=-zt1
            zt2=uelb(i)*gmatrix_3(i1,j,i2,i,l)*
     &           vpos(l)*(q2(j)-ppos(j))
            if(j.lt.4)zt2=-zt2

            zout1=zout1+zt1
            zout2=zout2+zt2

         enddo

         zout1=zout1+uelb(i)*mq*gmatrix_2(i2,i1,i,l)*vpos(l)
         zout2=zout2+uelb(i)*mq*gmatrix_2(i1,i2,i,l)*vpos(l)

 900  enddo


      zout1=zout1/(ul-mq**2)
      zout2=zout2/(tl-mq**2)

      zout=zout1+zout2
      zout=zout*4d0*pi*dsqrt(alphaEM(qsq1)*alphaEM(qsq2))
      zout=zout*dsqrt(conv)

      zoutarr(p,i1,i2)=zout

      enddo
      enddo

      return
      end
