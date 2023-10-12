ccc   calculates proton opacity
      subroutine opacity(i,j,bt,out,out1)
      implicit double precision(a-y)
      integer i,j,it,nt

      include 'nchan.f'
      include 'pi.f'
      include 'vars.f'
      include 'survpars.f'

      ampi=0.02d0
      amro=1d0
      a4=4d0*ampi
      alo=dlog(amro/ampi)
      bpol=2.4d0

      nt=6000
      htt=6d0/dble(nt)**2

         out=0d0
         out1=0d0

         do it=0,nt
            t=dble(it)**2*htt
            if(it.eq.0) t=1d-8
            wt=htt*2d0*dble(it)/4d0/pi
            if(it.eq.0) wt=wt/2d0
            bes0=besj0(bt*dsqrt(t))

            ffi=dexp(-((t+0.08d0+bb0(i))*bex(i))**cc0(i)+
     &           (bex(i)*(bb0(i)+0.08d0))**cc0(i))
            ffj=dexp(-((t+0.08d0+bb0(j))*bex(j))**cc0(j)+
     &           (bex(j)*(bb0(j)+0.08d0))**cc0(j))

            asp1=asp
            form1=dlog(ffi*ffj)-2d0*t*asp1*dlog(rts)
cccccccccccccccc
            ah=dsqrt(1d0+a4/t)
            h1pi=2d0*a4+t*(alo-ah*ah*ah*dlog((1d0+ah)/(ah-1d0)))
            h1pi=h1pi*sigo/(72d0*pi**3*(1d0+t/bpol)**2)
ccccccccccccccc
            ww=bes0*dexp(form1-2d0*h1pi*dlog(rts))
            aspt=t*asp+h1pi

            out=out+ww*wt
            out1=out1+bes0*wt*ffi*ffj

         enddo

      return
      end
