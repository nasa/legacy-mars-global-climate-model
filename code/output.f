      SUBROUTINE OUTPUT(TAU,ROTPER,RESTRT,IDAY,KEY8,KEY9)

!     Legacy Mars GCM v24
!     Mars Climate Modeling Center
!     NASA Ames Research Center

C  PURPOSE:
C      OUTPUT IS A SKELETAL ROUTINE THAT CAN BE ADAPTED FOR PRODUCING
C      OUTPUT OTHER THAN HISTORY TAPE OUTPUT.
C  MODIFIED FROM ORIGINAL BY
C      STEVE POHORSKY    INFORMATICS     TASK 605    JUL 82
C  FOR
C      JIM POLLACK
C  ENVIRONMENT
C      Cray-2            UNICOS 3.0      FORTRAN
C  REVISION HISTORY
C      Re-written September 1993.
C      Removed references to all commons and created an argument list
C      to pass all variables.
C
C  INPUT PARAMETERS
C      TAU      - CURRENT SIMULATED TIME MEASURED IN HOURS FROM WHEN THE
C                 COLD START (DAY 0 AND HOUR 0) WAS PERFORMED.
C      TOFDAY   - TIME IN (MARTIAN) HOURS MEASURED FROM THE START OF THE
C                 CURRENT SIMULATED DAY.
C      RESTRT   - TRUE UNTIL THE INITIALIZATION SECTION OF THIS PROGRAM
C                 FINISHES.
C      KEY8 AND - THESE FLAGS ARE TRUE ONCE EVERY TAUO AND TAUO2
C       KEY9      (RESPECTIVELY) HOURS. THESE CAN BE USED TO OUTPUT
C                 DIFFERENT QUANTITIES AT DIFFERENT FREQUENCIES.
C  CALLED BY
C      MAIN
C
      use grid_h
      use defines_h

      implicit none

C######################################################################

      LOGICAL KEY8, KEY9, RESTRT
      integer :: iday
      real*8  :: tau, rotper

C#=====================================================================

      IDAY = TAU/ROTPER

      IF(.NOT.RESTRT) THEN
        KEY8 = .FALSE.
        KEY9 = .FALSE.
      END IF

      RETURN

      END
