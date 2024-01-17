ccc   interpolator for Sudakov factor
      subroutine sudint(qtsq,mxx,tg,dtg)
      implicit double precision(a-y)
      integer i,j
      double precision iqinc,iminc

      include 'sudpars.f'

      mx=dlog(mxx)
      iqinc=1.0d0/qinc 
      iminc=1.0d0/minc 
      i=max(1,ceiling((qtsq-qmin)*iqinc))
      j=max(1,ceiling((mx-mmin)*iminc))


cccccccccccc

      m1=(tgint(3,i+1,j)-tgint(3,i,j))*iqinc
      tg1=tgint(3,i,j)+m1*(qtsq-tgint(1,i,j))

      m1=(tgint(4,i+1,j)-tgint(4,i,j))*iqinc
      dtg1=tgint(4,i,j)+m1*(qtsq-tgint(1,i,j))

      m2=(tgint(3,i+1,j+1)-tgint(3,i,j+1))*iqinc
      tg2=tgint(3,i,j+1)+m2*(qtsq-tgint(1,i,j+1))

      m2=(tgint(4,i+1,j+1)-tgint(4,i,j+1))*iqinc
      dtg2=tgint(4,i,j+1)+m2*(qtsq-tgint(1,i,j+1))

ccccccccccccc

      mf=(tg2-tg1)*iminc
      tg=tg1+mf*(mx-tgint(2,i,j))

      mf=(dtg2-dtg1)*iminc
      dtg=dtg1+mf*(mx-tgint(2,i,j))

      return
      end
