      subroutine wwpol_match(p,ip1,ip2,match)
      implicit none
      logical match
      integer p,ip1,ip2

      match=.false.
      
      if(p.eq.1)then
         if(ip1.eq.1.and.ip2.eq.1)match=.true.
      elseif(p.eq.2)then
         if(ip1.eq.1.and.ip2.eq.-1)match=.true.
      elseif(p.eq.3)then
         if(ip1.eq.-1.and.ip2.eq.1)match=.true.
      elseif(p.eq.4)then
         if(ip1.eq.-1.and.ip2.eq.-1)match=.true.
      elseif(p.eq.5)then
         if(ip1.eq.0.and.ip2.eq.1)match=.true.
      elseif(p.eq.6)then
         if(ip1.eq.0.and.ip2.eq.-1)match=.true.
      elseif(p.eq.7)then
         if(ip1.eq.1.and.ip2.eq.0)match=.true.
      elseif(p.eq.8)then
         if(ip1.eq.-1.and.ip2.eq.0)match=.true.
      elseif(p.eq.9)then
         if(ip1.eq.0.and.ip2.eq.0)match=.true.
      endif

c$$$      if(p.eq.1)then
c$$$         match=.true.
c$$$      else
c$$$         match=.false.
c$$$      endif

c      match=.true.
      
      return
      end
