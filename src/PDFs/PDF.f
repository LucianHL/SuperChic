ccc   calculates single PDF xg and dxg/dQ^2
      function xg(x,qsq)
      implicit double precision (a-z)
      double precision garr(-6:6)

      q0=1.5d0
      eps=1d-3
      eps1=1d-4

      xg=0d0
      
      if(qsq.gt.q0)then

         glu=xgi(x,qsq)

      else       
         
         qp=q0+eps1
         qm=q0-eps1

         glu1p=xgi(x,qp)
         if(glu1p.lt.1d-20)return
         glu1pp=xgi(x,qp+eps)
         glu1m=xgi(x,qm)
         glu1mp=xgi(x,qm+eps)
         
         diffp=qp*(glu1pp-glu1p)/(eps*glu1p)
         diffm=qm*(glu1mp-glu1m)/(eps*glu1m)

         ddiff=q0*(diffp-diffm)/(2d0*eps1)

         lamn=((diffp+2d0*dlog(q0))/(1d0+2d0*dlog(q0))-(ddiff+2d0*
     &        dlog(q0)+4d0)/4d0/(1d0+dlog(q0)))
         lamn=lamn/((1d0+dlog(q0))/(1d0+2d0*dlog(q0))-
     &        (2d0+dlog(q0))/4d0/(1d0+dlog(q0)))
         
         lamt1=(diffp+2d0*dlog(q0)-lamn*(1d0+dlog(q0)))
     &        /(1d0+2d0*dlog(q0))
         lamt2=(ddiff+2d0*dlog(q0)+4d0-lamn*(2d0+dlog(q0)))
     &        /4d0/(1d0+dlog(q0))
   
         power=2d0+(lamn-2d0)*qsq/(q0)+lamt2*qsq**2/q0**2

         if(power.lt.0d0)then
            glu=glu1p
         else
            call fglow1(q0,lamn,lamt2,qsq,glu1p,glu)
         endif

      endif
        
      xg=glu

      if(x.gt.0.9999d0)then
         xg=0d0
      endif

      return
      end

      function dxg(x,qsq)
      implicit double precision (a-z)
      double precision garr(-6:6)

      q0=1.5d0
      eps=1d-3
      eps1=1d-4

      dxg=0d0

      if(qsq.gt.q0)then

      glu=xgi(x,qsq)
      glup=xgi(x,qsq+eps)

      dxg=(glup-glu)/eps

      else

        qp=q0+eps1
         qm=q0-eps1

         glu1p=xgi(x,qp)
         if(glu1p.lt.1d-20)return
         glu1pp=xgi(x,qp+eps)
         glu1m=xgi(x,qm)
         glu1mp=xgi(x,qm+eps)

         diffp=qp*(glu1pp-glu1p)/(eps*glu1p)
         diffm=qm*(glu1mp-glu1m)/(eps*glu1m)
         ddiff=q0*(diffp-diffm)/(2d0*eps1)

         lamn=((diffp+2d0*dlog(q0))/(1d0+2d0*dlog(q0))-(ddiff+2d0*
     &        dlog(q0)+4d0)/4d0/(1d0+dlog(q0)))
         lamn=lamn/((1d0+dlog(q0))/(1d0+2d0*dlog(q0))-
     &        (2d0+dlog(q0))/4d0/(1d0+dlog(q0)))
         
         lamt1=(diffp+2d0*dlog(q0)-lamn*(1d0+dlog(q0)))
     &        /(1d0+2d0*dlog(q0))
         lamt2=(ddiff+2d0*dlog(q0)+4d0-lamn*(2d0+dlog(q0)))
     &        /4d0/(1d0+dlog(q0))
         
         power=2d0+(lamn-2d0)*qsq/(q0)+lamt2*qsq**2/q0**2

         if(power.lt.0d0)then
            dxg=0d0
         else
            call fglow1(q0,lamn,lamt2,qsq,glu1p,glu)
            dxg=((lamn-2d0)/q0+2d0*lamt2*qsq/q0**2)*dlog(qsq)
            dxg=dxg+power/qsq
            dxg=dxg*glu
         endif

      endif

      if(x.gt.0.9999d0)then
         dxg=0d0
      endif
         
      return
      end

      subroutine fglow1 (q0,lam,lamt,q,gluin,glout)
      implicit double precision (a-z)

      a=gluin/(q0**(lam+lamt))
      glout=a*q**(2d0+(lam-2d0)*q/(q0)+lamt*q**2/q0**2)

      return
      end

