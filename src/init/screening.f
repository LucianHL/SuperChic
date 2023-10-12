ccc   calculates screened amplitude (in kt space)
      subroutine screening(i,j,ktsq,out,out1)
      implicit double precision(a-y)
      integer i,j,ib,nb

      include 'nchan.f'
      include 'pi.f'
      include 'vars.f'
      include 'survpars.f'

      nb=5000
      hb=99d0/dble(nb)  

      out=0d0
      out1=0d0

      do ib=1,nb
         bt=dble(ib)*hb   
         wt=-bt/2d0/pi*hb
 
         call opacityint(i,j,bt,fr,fr1)
      
         sige=sigo*dexp(dlog(rts)*2d0*ep)
         fr=fr*gaa(i)*gaa(j)*sige
         fr1=fr1*gaa(i)*gaa(j)*sige
         
         out=out+wt*(1d0-dexp(-fr/2d0))*besj0(bt*dsqrt(ktsq))
     &        *gaa(i)*gaa(j)
         out1=out1+wt*(1d0-dexp(-fr1/2d0))*besj0(bt*dsqrt(ktsq))
     &        *gaa(i)*gaa(j)

      enddo

      return
      end
