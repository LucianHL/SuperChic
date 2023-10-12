ccc   returns length of character
      subroutine length(char,l1)
      character *(*) char
      integer i,j,l1

      do i=len(char),1,-1
      if(char(i:i).eq.' ')then
         j=i
      else
         goto 777
      endif
     
      enddo

 777  l1=j-1

      return
      end


