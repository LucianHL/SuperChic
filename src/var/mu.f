ccc   sets hard scale
      subroutine setmu(mu)
      implicit none
      double precision mu

      include 'mt.f'
      include 'vars.f'
      include 'mfact.f'

      if(mfact.eq.'mx')then
         mu=mx
      elseif(mfact.eq.'mt')then
         mu=mt*2d0
      endif

      return
      end
