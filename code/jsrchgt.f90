      function jsrchgt(N,SX,INC,TARGET)

!     Legacy Mars GCM v24
!     Mars Climate Modeling Center
!     NASA Ames Research Center

!      Find the first array element that is greater than the 
!      TARGET value.  If N < 1, then 0 is returned.  If no
!      value is found, N+1 is returned.  Replaces Cray version
!      on the workstation.
!      Bisection search implemented 09/29/93 to improve speed.
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
!                   greater than the TARGET value.
!----------------------------------------------------------------------C

      implicit none
 
      integer :: n, inc, jsrchgt, ians, jl, jh, jm
      real*8  :: SX(N), target

!======================================================================C

      if(N.lt.1) then
        ians = 0
      elseif(TARGET.gt.SX(N)) then
        ians = N+1
      elseif(TARGET.lt.SX(1)) then
        ians = 1
      else

        JL = 1
        JH = N

   10   CONTINUE
        if(JH-JL.gt.1) then
          JM = (JL+JH)/2

          if(TARGET.GT.SX(JM)) then
            JL = JM
            JM = (JL+JH)/2
          else
            JH = JM
            JM = (JL+JH)/2
          end if

          GOTO 10
        end if

        if(TARGET.EQ.SX(JH)) JH = JH+1

        ians = JH
      end if

      JSRCHGT = ians

      end function jsrchgt
