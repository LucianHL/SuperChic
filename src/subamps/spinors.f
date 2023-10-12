      subroutine up(i,out)
      implicit none
      double precision phiu,phi,thetau,pmod,norm,cthetau
      complex*16 out(4)
      integer i
      
      include 'mom.f'
      include 'zi.f'
      include 'mq.f'

      pmod=dsqrt(q(1,i)**2+q(2,i)**2+q(3,i)**2)
      norm=dsqrt(q(4,i)+mq)
      
      cthetau=q(3,i)/pmod
      thetau=dacos(cthetau)
      phiu=phi(i)

      out(1)=norm*dcos(thetau/2d0)
      out(2)=norm*dsin(thetau/2d0)*(dcos(phiu)+zi*dsin(phiu))
      out(3)=pmod*dcos(thetau/2d0)/norm
      out(4)=pmod*dsin(thetau/2d0)/norm*(dcos(phiu)+zi*dsin(phiu))
      
      return
      end

      subroutine um(i,out)
      implicit none
      double precision phiu,phi,thetau,pmod,norm,cthetau
      complex*16 out(4)
      integer i
      
      include 'mom.f'
      include 'zi.f'
      include 'mq.f'

      pmod=dsqrt(q(1,i)**2+q(2,i)**2+q(3,i)**2)
      norm=dsqrt(q(4,i)+mq)
      
      cthetau=q(3,i)/pmod
      thetau=dacos(cthetau)
      phiu=phi(i)
      
      out(1)=-norm*dsin(thetau/2d0)
      out(2)=norm*dcos(thetau/2d0)*(dcos(phiu)+zi*dsin(phiu))
      out(3)=pmod*dsin(thetau/2d0)/norm
      out(4)=-pmod*dcos(thetau/2d0)/norm*(dcos(phiu)+zi*dsin(phiu))
      
      return
      end

      subroutine vp(i,out)
      implicit none
      double precision phiu,phi,thetau,pmod,norm,cthetau
      complex*16 out(4)
      integer i
      
      include 'mom.f'
      include 'zi.f'
      include 'mq.f'

      pmod=dsqrt(q(1,i)**2+q(2,i)**2+q(3,i)**2)
      norm=dsqrt(q(4,i)+mq)
      
      cthetau=q(3,i)/pmod
      thetau=dacos(cthetau)
      phiu=phi(i)
      
      out(1)=pmod*dsin(thetau/2d0)/norm
      out(2)=-pmod*dcos(thetau/2d0)*(dcos(phiu)+zi*dsin(phiu))/norm
      out(3)=-norm*dsin(thetau/2d0)
      out(4)=norm*dcos(thetau/2d0)*(dcos(phiu)+zi*dsin(phiu))
      
      return
      end

      subroutine vm(i,out)
      implicit none
      double precision phiu,phi,thetau,pmod,norm,cthetau
      complex*16 out(4)
      integer i
      
      include 'mom.f'
      include 'zi.f'
      include 'mq.f'
      
      pmod=dsqrt(q(1,i)**2+q(2,i)**2+q(3,i)**2)
      norm=dsqrt(q(4,i)+mq)
      
      cthetau=q(3,i)/pmod
      thetau=dacos(cthetau)
      phiu=phi(i)

      out(1)=pmod*dcos(thetau/2d0)/norm
      out(2)=pmod*dsin(thetau/2d0)*(dcos(phiu)+zi*dsin(phiu))/norm
      out(3)=norm*dcos(thetau/2d0)
      out(4)=norm*dsin(thetau/2d0)*(dcos(phiu)+zi*dsin(phiu))
      
      return
      end

      subroutine upb(i,outb)
      implicit none
      complex*16 outb(4),out(4)
      integer i,j,k
      
      include 'gmatrices.f'

      call up(i,out)

      do j=1,4
         outb(j)=0d0
      enddo
      
      do j=1,4
         do k=1,4
            outb(j)=outb(j)+gmatrix(4,j,k)*dconjg(out(k))
         enddo
      enddo

      return
      end

      subroutine umb(i,outb)
      implicit none
      complex*16 outb(4),out(4)
      integer i,j,k
      
      include 'gmatrices.f'

      call um(i,out)

      do j=1,4
         outb(j)=0d0
      enddo
      
      do j=1,4
         do k=1,4
            outb(j)=outb(j)+gmatrix(4,j,k)*dconjg(out(k))
         enddo
      enddo

      return
      end
     
      subroutine vpb(i,outb)
      implicit none
      complex*16 outb(4),out(4)
      integer i,j,k
      
      include 'gmatrices.f'

      call vp(i,out)

      do j=1,4
         outb(j)=0d0
      enddo
      
      do j=1,4
         do k=1,4
            outb(j)=outb(j)+gmatrix(4,j,k)*dconjg(out(k))
         enddo
      enddo

      return
      end

      subroutine vmb(i,outb)
      implicit none
      complex*16 outb(4),out(4)
      integer i,j,k
      
      include 'gmatrices.f'

      call vm(i,out)

      do j=1,4
         outb(j)=0d0
      enddo
      
      do j=1,4
         do k=1,4
            outb(j)=outb(j)+gmatrix(4,j,k)*dconjg(out(k))
         enddo
      enddo

      return
      end

      subroutine slash(p,out)
      implicit none
      double precision p(4)
      complex*16 out(4,4)
      integer i,j,k
      
      include 'gmatrices.f'

      do i=1,4
         do j=1,4
            out(i,j)=0d0
         enddo
      enddo
            
      do k=1,4
         do i=1,4
            do j=1,4
               if(k.lt.4)then
                  out(i,j)=out(i,j)-p(k)*gmatrix(k,i,j)
               else
                  out(i,j)=out(i,j)+p(k)*gmatrix(k,i,j)
               endif
            enddo
         enddo
      enddo
               
      return
      end
