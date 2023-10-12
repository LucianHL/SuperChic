ccc   writes proton opacity to array 
      subroutine calcop
      implicit double precision(a-y)
      integer nb,i1,i2,ib

      include 'nchan.f'
      include 'opac.f'

      nb=900
      hb=100d0/dble(nb)

      print*,'Calculating opacity...'

      do ib=1,nb+1
         bt=dble(ib-1)*hb

         do i1=1,nch
            do i2=1,nch

               call opacity(i1,i2,bt,fr,fr1)

               op(i1,i2,ib,1)=bt
               op(i1,i2,ib,2)=fr
               oph(i1,i2,ib,1)=bt
               oph(i1,i2,ib,2)=fr1

            enddo
         enddo
      enddo

      print*,'Done!'

      return
      end
