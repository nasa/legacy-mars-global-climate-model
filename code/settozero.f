       subroutine settozero(nsize,array)

!     Legacy Mars GCM v24
!     Mars Climate Modeling Center
!     NASA Ames Research Center

       implicit none

       integer nsize
       real*8 array(nsize)

       integer i

       do i=1,nsize
         array(i) = 0.
       enddo

       return

       end
       subroutine settozero4(nsize,array)

       implicit none

       integer nsize
       real*4 array(nsize)

       integer i

       do i=1,nsize
         array(i) = 0.0
       enddo

       return

       end
