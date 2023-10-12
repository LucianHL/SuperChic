ccc   interpolator for skewed pdf
      subroutine hpdfint(x,qtsq,tg,dtg)
      implicit double precision(a-y)
      integer i,j

      include 'hgpars.f'

      i=nint((qtsq-qmin)/qinc)
      j=nint((x-lxmin)/lxinc)

      if(dble(i).lt.(qtsq-qmin)/qinc)then
         i=i+1
      endif

      if(dble(j).lt.(x-lxmin)/lxinc)then
         j=j+1
      endif

      if(i.eq.0)i=i+1
      if(j.eq.0)j=j+1

cccccccccccc  
      
      m1=(hgint(3,i+1,j)-hgint(3,i,j))/qinc
      tg1=hgint(3,i,j)+m1*(qtsq-hgint(1,i,j))

      tg1p=hgint(3,i+1,j)+m1*(qtsq-hgint(1,i+1,j))

      m1=(hgint(4,i+1,j)-hgint(4,i,j))/qinc
      dtg1=hgint(4,i,j)+m1*(qtsq-hgint(1,i,j))

      m2=(hgint(3,i+1,j+1)-hgint(3,i,j+1))/qinc
      tg2=hgint(3,i,j+1)+m2*(qtsq-hgint(1,i,j+1))

      m2=(hgint(4,i+1,j+1)-hgint(4,i,j+1))/qinc
      dtg2=hgint(4,i,j+1)+m2*(qtsq-hgint(1,i,j+1))
      
ccccccccccccc

      mf=(tg2-tg1)/lxinc
      tg=tg1+mf*(x-hgint(2,i,j))
    
      mf=(dtg2-dtg1)/lxinc
      dtg=dtg1+mf*(x-hgint(2,i,j))

      return
      end
