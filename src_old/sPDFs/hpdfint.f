ccc   interpolator for skewed pdf
      subroutine hpdfint(x,qtsq,tg,dtg)
      implicit double precision(a-y)
      integer i,j
      double precision iqinc, ilxinc

      include 'hgpars.f'
      iqinc=1.0d0/qinc 
      ilxinc=1.0d0/lxinc 
       
      i=max(1,ceiling((qtsq-qmin)*iqinc))
      j=max(1,ceiling((x-lxmin)*ilxinc))


cccccccccccc

      m1=(hgint(3,i+1,j)-hgint(3,i,j))*iqinc
      tg1=hgint(3,i,j)+m1*(qtsq-hgint(1,i,j))

      tg1p=hgint(3,i+1,j)+m1*(qtsq-hgint(1,i+1,j))

      m1=(hgint(4,i+1,j)-hgint(4,i,j))*iqinc
      dtg1=hgint(4,i,j)+m1*(qtsq-hgint(1,i,j))

      m2=(hgint(3,i+1,j+1)-hgint(3,i,j+1))*iqinc
      tg2=hgint(3,i,j+1)+m2*(qtsq-hgint(1,i,j+1))

      m2=(hgint(4,i+1,j+1)-hgint(4,i,j+1))*iqinc
      dtg2=hgint(4,i,j+1)+m2*(qtsq-hgint(1,i,j+1))

ccccccccccccc

      mf=(tg2-tg1)*ilxinc
      tg=tg1+mf*(x-hgint(2,i,j))

      mf=(dtg2-dtg1)*ilxinc
      dtg=dtg1+mf*(x-hgint(2,i,j))

      return
      end
