c************************************************************************
c**                                                                    **
c**           MadGraph/MadEvent Interface to FeynRules                 **
c**                                                                    **
c**          C. Duhr (Louvain U.) - M. Herquet (NIKHEF)                **
c**                                                                    **
c************************************************************************

      subroutine setpara(param_name)
      implicit none

      character*(*) param_name
      logical readlha

      include 'coupl.inc'
      include 'input.inc'
      include 'model_functions.inc'

      integer maxpara
      parameter (maxpara=5000)
      
      integer npara
      character*20 param(maxpara),value(maxpara)

      logical updateloop
      common /to_updateloop/updateloop
      data updateloop /.true./

      call LHA_loadcard(param_name,npara,param,value)
      ! also loop parameters should be initialised here
      if (updateloop) then
         include 'param_read.inc'
         call coup()         
      else   
         updateloop=.true.
         include 'param_read.inc'
         call coup()
         updateloop=.false.
      endif
      return

      end
c$$$
c$$$      subroutine setParamLog(OnOff)
c$$$
c$$$      logical OnOff
c$$$      logical WriteParamLog
c$$$      data WriteParamLog/.TRUE./
c$$$      common/IOcontrol/WriteParamLog
c$$$
c$$$      WriteParamLog = OnOff
c$$$
c$$$      end
c$$$
c$$$      subroutine setpara2(param_name)
c$$$      implicit none
c$$$
c$$$      character(512) param_name
c$$$
c$$$      integer k
c$$$      logical found
c$$$
c$$$      character(512) ParamCardPath
c$$$      common/ParamCardPath/ParamCardPath
c$$$
c$$$      if (param_name(1:1).ne.' ') then
c$$$        ! Save the basename of the param_card for the ident_card.
c$$$        ! If no absolute path was used then this ParamCardPath
c$$$        ! remains empty
c$$$        ParamCardPath = '.'
c$$$        k = LEN(param_name)
c$$$        found = .False.
c$$$        do while (k.ge.1.and..not.found)
c$$$          if (param_name(k:k).eq.'/') then
c$$$              found=.True.
c$$$          endif
c$$$          k=k-1
c$$$        enddo
c$$$        if (k.ge.1) then
c$$$          ParamCardPath(1:k)=param_name(1:k)
c$$$        endif
c$$$        call setpara(param_name)
c$$$      endif
c$$$      if (param_name(1:1).eq.'*') then
c$$$         ! Dummy call to printout so that it is available in the
c$$$         ! dynamic library for MadLoop BLHA2
c$$$         ! In principle the --whole-archive option of ld could be 
c$$$         ! used but it is not always supported
c$$$         call printout()
c$$$         call setParamLog(.True.)
c$$$      endif
c$$$      return
c$$$
c$$$      end
