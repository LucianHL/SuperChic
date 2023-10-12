ccc   integrates gg-->MM amplitude over momentum fractions x,y
ccc   and writes to array
      subroutine calcmes
      implicit double precision(a-y)
      integer i,itot,jx,jxtot,jy,jytot,p

      include 'mes.f'
      include 'polarization.f'
      include 'proc.f'
      include 'mixing.f'

      print*,'Integrating meson pair amplitude over x,y...'

      jxtot=1000
      jytot=1000
      jxinc=1d0/dble(jxtot)
      jyinc=1d0/dble(jytot)
      
      cmax=0.99d0
      cmin=0d0
      cmin=-0.99d0

      itot=200
      cinc=(cmax-cmin)/dble(itot)

      do p=1,pol

      do i=1,itot+1
         
         cost=cmin+cinc*dble(i-1)

         sum=0d0
         sum1=0d0
         sum2=0d0
         sum3=0d0
         sum4=0d0
         sum5=0d0

         do 555 jx=1,jxtot
         do 555 jy=1,jytot
            
            x=(dble(jx)-0.5d0)*jxinc
            y=(dble(jy)-0.5d0)*jyinc

            if(proc.eq.9.or.proc.eq.10.or.proc.eq.11.or.proc.eq.12)then
               call pipixy(x,y,cost,out)
            elseif(proc.eq.13.or.proc.eq.17)then
               call rhorhoxy(p,x,y,cost,out)
            elseif(proc.eq.14.or.proc.eq.15.or.proc.eq.16)then
               call pipixy(x,y,cost,out)
            endif

            g2x=(-1.5d0+7.5d0*(2d0*x-1d0)**2)
            g4x=15d0/8d0*(1d0-14d0*x**2+21d0*x**4)
            g2y=(-1.5d0+7.5d0*(2d0*y-1d0)**2)
            g4y=15d0/8d0*(1d0-14d0*y**2+21d0*y**4)

            sum=sum+out*jxinc*jyinc
            sum1=sum1+out*jxinc*jyinc*a28*g2x
            sum2=sum2+out*jxinc*jyinc*a28**2*g2x*g2y

            sum3=sum3+out*jxinc*jyinc*a48*g4x
            sum4=sum4+out*jxinc*jyinc*a28*a48*g2y*g4x
            sum5=sum5+out*jxinc*jyinc*a48**2*g4x*g4y

 555     enddo

         mesamp(p,1,1,i)=cost
         mesamp(p,1,2,i)=cost
         mesamp(p,1,3,i)=cost
         mesamp(p,1,4,i)=cost
         mesamp(p,1,5,i)=cost
         mesamp(p,1,6,i)=cost
   
         mesamp(p,2,1,i)=sum
         mesamp(p,2,2,i)=sum1
         mesamp(p,2,3,i)=sum2
         mesamp(p,2,4,i)=sum3
         mesamp(p,2,5,i)=sum4
         mesamp(p,2,6,i)=sum5
     
      enddo

      enddo

      print*,'Done!'

      return
      end
