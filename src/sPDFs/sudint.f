ccc   interpolator for Sudakov factor
      subroutine sudint(qtsq,mxx,tg,dtg)
      implicit double precision(a-y)
      integer i,j

      include 'sudpars.f'

      mx=dlog(mxx)

      i=nint((qtsq-qmin)/qinc)
      j=nint((mx-mmin)/minc)

      if(dble(i).lt.(qtsq-qmin)/qinc)then
         i=i+1
      endif

      if(dble(j).lt.(mx-mmin)/minc)then
         j=j+1
      endif

      if(i.eq.0)i=i+1
      if(j.eq.0)j=j+1

cccccccccccc  
      
      m1=(tgint(3,i+1,j)-tgint(3,i,j))/qinc
      tg1=tgint(3,i,j)+m1*(qtsq-tgint(1,i,j))

      m1=(tgint(4,i+1,j)-tgint(4,i,j))/qinc
      dtg1=tgint(4,i,j)+m1*(qtsq-tgint(1,i,j))

      m2=(tgint(3,i+1,j+1)-tgint(3,i,j+1))/qinc
      tg2=tgint(3,i,j+1)+m2*(qtsq-tgint(1,i,j+1))

      m2=(tgint(4,i+1,j+1)-tgint(4,i,j+1))/qinc
      dtg2=tgint(4,i,j+1)+m2*(qtsq-tgint(1,i,j+1))
      
ccccccccccccc

      mf=(tg2-tg1)/minc
      tg=tg1+mf*(mx-tgint(2,i,j))

      mf=(dtg2-dtg1)/minc
      dtg=dtg1+mf*(mx-tgint(2,i,j))

      return
      end
