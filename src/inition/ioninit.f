      subroutine ioninit
      implicit none

      include 'pi.f'
      include 'vars.f'
      include 'intag.f'
      include 'pdfinf.f'
      include 'ionqcd.f'
      include 'beam.f'
      include 'qcd.f'
      include 'ion.f'
      include 'p0Xn.f'

      
      pi=dacos(-1d0)

cccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccc

      print*,'Initialising ion screening effects...'

      if(beam.eq.'ionp'.and.qcd.eqv..true.)ionqcd='coh'
      if(faa.eq.'0X'.or.faa.eq.'X0')then
         if(ionbreakup)pAAvar=.true.
      endif
      
      call ionpars
      

      if(ionbreakup)call gdrin
      call opacpars(rtsnn)
      call rhonorm      
      call rhoxycalc      
      call tpcalc
      print*,'gdrset...'
      if(ionbreakup)call gdrset      
      print*,'...done'

      call opacpcalc
      if(beam.eq.'ion')then

         if(pAAvar)then
         
            do ifaa=1,3
c               print*,'ifaa = ',ifaa
               call opacpbcalc
               call screencalc
            enddo
         
         return
         endif


         print*,'opapbcalc...'
         call opacpbcalc
         print*,'...done'
         if(ionqcd.eq.'incoh')then
            call tpqcdcalc
            call s2qcdion
         endif
      elseif(beam.eq.'ionp')then
         call opacpbpcalc
         if(ionqcd.eq.'incoh')then
            call tpqcdcalc
            call s2qcdionp
         endif
      endif
      if(qcd)then
         if(ionqcd.eq.'coh')call screencalc
      else
         call screencalc
      endif
      
      print*,'...finished!'

      return
      end
