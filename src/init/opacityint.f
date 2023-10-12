ccc   interpolates opacity from array
      subroutine opacityint(i,j,bt,out,out1)
      implicit double precision(a-y)
      integer i,j,it

      include 'opac.f'

      incbt=op(1,1,2,1)-op(1,1,1,1)
      it=nint(bt/incbt)
      if(dble(it).gt.(bt/incbt))then
         it=it-1
      endif

      m=(op(i,j,it+2,2)-op(i,j,it+1,2))
     &/(op(i,j,it+2,1)-op(i,j,it+1,1))
      del=bt-op(1,1,it+1,1)
      mh=(oph(i,j,it+2,2)-oph(i,j,it+1,2))
     &/(oph(i,j,it+2,1)-oph(i,j,it+1,1))
      delh=bt-oph(1,1,it+1,1)

      out=m*del+op(i,j,it+1,2)
      out1=mh*delh+oph(i,j,it+1,2)

      return
      end
