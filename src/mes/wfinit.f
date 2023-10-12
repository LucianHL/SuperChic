ccc   initialises parameters for meson wave function evolution
      subroutine wfinit
      implicit double precision(a-y)
      integer i,j

      include 'anom.f'

      cf=4d0/3d0
      nf=3d0
      nc=3d0
      beta0=11d0-2d0*nf/3d0

      do i=2,6,2

         ii=dble(i)

         num1=0d0
         do j=1,i+1
            num1=num1+1d0/dble(j)
         enddo

         gamqq(i)=cf*(3d0+2d0/(ii+1d0)/(ii+2d0)-4d0*num1)
         gamgg(i)=beta0+nc*(8d0/(ii+1d0)/(ii+2d0)-4d0*num1)
         gamgq(i)=nf*12d0/(ii+1d0)/(ii+2d0)
         gamqg(i)=cf*ii*(ii+3d0)/3d0/(ii+1d0)/(ii+2d0)

         gamp(i)=(gamqq(i)+gamgg(i)+dsqrt((gamqq(i)-gamgg(i))**2+
     &        4d0*gamqg(i)*gamgq(i)))/2d0
         gamm(i)=(gamqq(i)+gamgg(i)-dsqrt((gamqq(i)-gamgg(i))**2+
     &        4d0*gamqg(i)*gamgq(i)))/2d0

         rhop(i)=6d0*gamgq(i)/(gamp(i)-gamgg(i))
         rhom(i)=gamqg(i)/(gamm(i)-gamqq(i))/6d0

      enddo

      return
      end
