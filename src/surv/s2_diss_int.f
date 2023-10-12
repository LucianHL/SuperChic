ccc   interpolator for Sudakov factor
      subroutine s2_diss_int(ii,qt,x,tg)
      implicit none
      double precision tg1,tg2,m1,m2,mf
      double precision lxinc,lxmin,lx,lqtmin,lqtinc,lqt
      double precision qt,x,tg
      integer i,j,ii
      integer nx,nqt

      include 's2_diss.f'

      if(ii.eq.3)then
         tg=s2dd        
         return
      endif
      
      nqt=10
      
      lx=dlog(x)
      lqt=dlog(qt)

      if(ii.eq.1)lqtmin=scearr(2,1,1)
      if(ii.eq.2)lqtmin=siearr(2,1,1)
      if(ii.eq.1)lxmin=scearr(1,1,1)
      if(ii.eq.2)lxmin=siearr(1,1,1)
      if(ii.eq.1)lqtinc=scearr(2,1,2)-scearr(2,1,1)
      if(ii.eq.2)lqtinc=siearr(2,1,2)-siearr(2,1,1)
      if(ii.eq.1)lxinc=scearr(1,2,1)-scearr(1,1,1)
      if(ii.eq.2)lxinc=siearr(1,2,1)-siearr(1,1,1)
      
      j=nint((lqt-lqtmin)/lqtinc)
      i=nint((lx-lxmin)/lxinc)

      if(dble(j).lt.(lqt-lqtmin)/lqtinc)then
         j=j+1
      endif

      if(dble(i).lt.(lx-lxmin)/lxinc)then
         i=i+1
      endif

      if(i.lt.1)i=1
      if(j.lt.1)j=1

c      print*,i,j

cccccccccccc  

      if(ii.eq.1)then

         if(j.gt.nqt)then
            
            m1=(scearr(3,i+1,nqt)-scearr(3,i,nqt))/lxinc
            tg=scearr(3,i,nqt)+m1*(lx-scearr(1,i,nqt))

         else
            
            m1=(scearr(3,i+1,j)-scearr(3,i,j))/lxinc
            tg1=scearr(3,i,j)+m1*(lx-scearr(1,i,j))
            
            m2=(scearr(3,i+1,j+1)-scearr(3,i,j+1))/lxinc
            tg2=scearr(3,i,j+1)+m2*(lx-scearr(1,i,j+1))
            
            mf=(tg2-tg1)/lqtinc
            tg=tg1+mf*(lqt-scearr(2,i,j))

         endif
            
      else

         if(j.gt.nqt)then
            
            m1=(siearr(3,i+1,nqt)-siearr(3,i,nqt))/lxinc
            tg=siearr(3,i,nqt)+m1*(lx-siearr(1,i,nqt))

         else
            
            m1=(siearr(3,i+1,j)-siearr(3,i,j))/lxinc
            tg1=siearr(3,i,j)+m1*(lx-siearr(1,i,j))
            
            m2=(siearr(3,i+1,j+1)-siearr(3,i,j+1))/lxinc
            tg2=siearr(3,i,j+1)+m2*(lx-siearr(1,i,j+1))
            
            mf=(tg2-tg1)/lqtinc
            tg=tg1+mf*(lqt-siearr(2,i,j))

         endif
            
      endif

      if(tg.lt.0d0)tg=0d0
      if(tg.gt.1d0)tg=1d0
      
      return
      end


      subroutine s2exp(qt,s2l,s2h,s2out)
      implicit none
      double precision s2out,qt,s2l,s2h,qt0,s2interp

      qt0=1d0
      s2out=s2l*s2interp(qt,qt0)+s2h*(1d0-s2interp(qt,qt0))

      return
      end

      
      subroutine s2exp2(qt1,qt2,s2ll,s2lh,s2hl,s2hh,s2out)
      implicit none
      double precision qt0,qt1,qt2,s2ll,s2lh,s2hl,s2hh,s2out
      double precision s2interp

      qt0=1d0
      s2out=s2ll*s2interp(qt1,qt0)*s2interp(qt2,qt0)
      s2out=s2out+s2hl*(1d0-s2interp(qt1,qt0))*s2interp(qt2,qt0)
      s2out=s2out+s2lh*(1d0-s2interp(qt2,qt0))*s2interp(qt1,qt0)
      s2out=s2out+s2hh*(1d0-s2interp(qt1,qt0))*(1d0-s2interp(qt2,qt0))

      return
      end

      function s2interp(qt,qt0)
      implicit none
      double precision s2interp,qt,qt0

      s2interp=dexp(-qt**2/qt0**2)

      return
      end
