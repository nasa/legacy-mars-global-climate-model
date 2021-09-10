       function jsrchge(N,SX,INC,TARGET)

!     Legacy Mars GCM v24
!     Mars Climate Modeling Center
!     NASA Ames Research Center

!      Find the first array element that is greater than or equal to
!      the TARGET value.  If N < 1, then 0 is returned.  If no value
!      is found, N+1 is returned.  Replaces Cray version on the 
!      workstation.
!
!      Started: 08/23/93
!
!      Input:
!        N      - Number of elements in array SX.
!        SX     - Array of numbers to be searched.  Assumed to be an
!                 ordered array.
!        INC    - Increment between elements of searched array.  Kept
!                 for compatibility with Cray call.
!        TARGET - Value searched for in array SX.
!
!      Output:
!        JSRCHGT - location in array SX where value of SX is first
!                   greater than or equal to the TARGET value.
!----------------------------------------------------------------------C

       implicit none

       integer :: n, inc, ians, jsrchge, i
       real*8  :: SX(N), target

!======================================================================C

       if(N.lt.1) then
         ians = 0
       else
         ians = N+1
         do I=1,N
           if(SX(I).GE.TARGET) then
             ians = I
             goto 10
           end if
         end do
       end if

   10 continue

      JSRCHGE = ians

      end function jsrchge

 
