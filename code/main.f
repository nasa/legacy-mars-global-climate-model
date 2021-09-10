C  Notices:

C  Copyright © 2021 United States Government as represented by the Administrator 
C  of the National Aeronautics and Space Administration.  All Rights Reserved.

C  Disclaimers

C  No Warranty: THE SUBJECT SOFTWARE IS PROVIDED "AS IS" WITHOUT ANY WARRANTY 
C  OF ANY KIND, EITHER EXPRESSED, IMPLIED, OR STATUTORY, INCLUDING, BUT NOT 
C  LIMITED TO, ANY WARRANTY THAT THE SUBJECT SOFTWARE WILL CONFORM TO 
C  SPECIFICATIONS, ANY IMPLIED WARRANTIES OF MERCHANTABILITY, FITNESS FOR A 
C  PARTICULAR PURPOSE, OR FREEDOM FROM INFRINGEMENT, ANY WARRANTY THAT THE SUBJECT 
C  SOFTWARE WILL BE ERROR FREE, OR ANY WARRANTY THAT DOCUMENTATION, IF PROVIDED, 
C  WILL CONFORM TO THE SUBJECT SOFTWARE. THIS AGREEMENT DOES NOT, IN ANY MANNER, 
C  CONSTITUTE AN ENDORSEMENT BY GOVERNMENT AGENCY OR ANY PRIOR RECIPIENT OF ANY RESULTS, 
C  RESULTING DESIGNS, HARDWARE, SOFTWARE PRODUCTS OR ANY OTHER APPLICATIONS RESULTING 
C  FROM USE OF THE SUBJECT SOFTWARE.  FURTHER, GOVERNMENT AGENCY DISCLAIMS ALL WARRANTIES 
C  AND LIABILITIES REGARDING THIRD-PARTY SOFTWARE, IF PRESENT IN THE ORIGINAL SOFTWARE, 
C  AND DISTRIBUTES IT "AS IS."

C  Waiver and Indemnity:  RECIPIENT AGREES TO WAIVE ANY AND ALL CLAIMS AGAINST THE 
C  UNITED STATES GOVERNMENT, ITS CONTRACTORS AND SUBCONTRACTORS, AS WELL AS ANY PRIOR 
C  RECIPIENT.  IF RECIPIENT'S USE OF THE SUBJECT SOFTWARE RESULTS IN ANY LIABILITIES, 
C  DEMANDS, DAMAGES, EXPENSES OR LOSSES ARISING FROM SUCH USE, INCLUDING ANY DAMAGES 
C  FROM PRODUCTS BASED ON, OR RESULTING FROM, RECIPIENT'S USE OF THE SUBJECT SOFTWARE, 
C  RECIPIENT SHALL INDEMNIFY AND HOLD HARMLESS THE UNITED STATES GOVERNMENT, ITS CONTRACTORS 
C  AND SUBCONTRACTORS, AS WELL AS ANY PRIOR RECIPIENT, TO THE EXTENT PERMITTED BY LAW.  
C  RECIPIENT'S SOLE REMEDY FOR ANY SUCH MATTER SHALL BE THE IMMEDIATE, UNILATERAL 
C  TERMINATION OF THIS AGREEMENT.

      PROGRAM MARSGCM

!     Legacy Mars GCM v24
!     Mars Climate Modeling Center
!     NASA Ames Research Center

C
C  PURPOSE:
C        THIS PROGRAM IS A NUMERICAL MODEL OF THE ATMOSPHERE OF MARS.
C        IT IS DESIGNED TO NUMERICALLY SIMULATE THE MARTIAN ATMOSPHERE
C        DURING A PERIOD OF TIME, TYPICALLY SEVERAL DAYS. IN THIS
C        PROGRAM THE ATMOSPHERE OF MARS IS REPRESENTED NUMERICALLY BY
C        A SET OF POINTS THAT FORMS A THREE-DIMENSIONAL GRID AROUND THE
C        SURFACE OF MARS. THE MODEL CALCULATES AND
C        OUTPUTS ATMOSPHERIC VARIABLES (WIND VELOCITY,
C        PRESSURE, TEMPERATURE, AND OTHER VALUES FOR EACH POINT OF THE
C        GRID) AT REGULAR TIME INTERVALS FOR THE PERIOD OF SIMULATED
C        TIME.
C
C        THIS PROGRAM HAS TWO BASIC PARTS, BOTH OF WHICH ARE CONTROLLED
C        BY THE MAIN PROGRAM MODULE.
C         1) MODEL INITIALIZATION
C             THIS PART IS EXECUTED FIRST AND IS DONE ONLY ONCE. ITS
C             PURPOSE IS TO SET UP FOR THE SECOND PART. MODEL
C             INITIALIZATION CONSISTS MAINLY OF READING INPUT VALUES,
C             CALCULATING CONSTANTS FOR USE IN THE SECOND PART, AND
C             INITIALIZING THE SET OF VALUES MAINTAINED FOR EACH GRID
C             POINT. THE INPUT VALUES DETERMINE WHETHER THE SETS OF
C             VALUES FOR THE GRID ARE INITIALIZED FROM THE OUTPUT FROM
C             A PREVIOUS RUN OF THIS PROGRAM OR FROM A SOMEWHAT
C             ARTIFICIAL SET OF VALUES (COLD START) CONTAINED IN THE
C             PROGRAM. THE MAJOR PROGRAM MODULES IN THIS PART ARE
C              INPUT, EMISS, MAGFAC, AND INSDET. INIT1 IS
C              ALSO USED IN THE CASE OF A COLD START.
C
C         2) CALCULATING ATMOSPHERIC VARIABLES AT SUCCESSIVE TIME
C            INTERVALS (TIME INTEGRATION)
C             THIS PART OF THE MODEL SOLVES, BY FINITE DIFFERENCE
C             TECHNIQUES, A NUMBER OF DIFFERENTIAL EQUATIONS THAT
C             EXPRESS THE PHYSICAL LAWS WHICH RELATE QUANTITIES SUCH
C             AS GROUND AND ATMOSPHERIC TEMPERATURES, PRESSURE, WIND
C             VELOCITY, AND THE AMOUNT OF HEATING BY THE SUN. THE
C             APPROACH USED IS TO INTEGRATE IN TIME, STARTING FROM THE
C             INITIAL VALUES SET BY THE FIRST PART, TO COMPUTE AT
C             REGULAR TIME INTERVALS THE SET OF VALUES MAINTAINED FOR
C             THE GRID. IN THE PROGRAM, THIS PART OF THE MODEL IS
C             PERFORMED AS A LOOP THAT RUNS FOR THE NUMBER OF TIME
C             INTERVALS (CALLED TIME STEPS) THAT ARE SPECIFIED FOR A
C             RUN THROUGH THE INPUT PARAMETERS. THE MAJOR PROGRAM
C             MODULES IN THIS PART ARE:  STEP, COMP1, COMP2, COMP3,
C             GMP, SDET, AND HISTORY.
C
C        FOR ADDITIONAL INFORMATION PLEASE SEE THE DOCUMENT 'THE MARS
C        GENERAL CIRCULATION MODEL (GCM) PROGRAM DESCRIPTION' AND
C        THE HEADER COMMENTS FOR THE OTHER PROGRAM MODULES.
C
C  MODIFIED FROM ORIGINAL BY
C      STEVE POHORSKY    INFORMATICS     TASK 605    SEP 82
C
C  FOR
C      JIM POLLACK
C  ENVIRONMENT
C      Cray-2           UNICOS 3.0         FORTRAN
C  REVISION HISTORY
C      March 1988    The entire GCM was converted from RATFOR to
C                    FORTRAN by JRS (Sterling Software Task 904)
C                    In this conversion NO attempt was made to fix
C                    up bad (read out-of-date) constructs or
C                    similar "improvements".  This was a straight
C                    change-this-RATFOR-do-loop-with ['s to a
C                    FORTRAN DO 789 type statement.  At a future
C                    date all routines should be redone to turn
C                    the code into a more uniformely STRUCTURED
C                    piece of work.
C
C  PROGRAM INPUT
C      TAPE4       - DATA FILE CONTAINING THE VALUES FOR THE TOPOG AND
C                    ALSP ARRAYS (TOPOGRAPHY AND SURFACE ALBEDO) IN
C                    UNFORMATTED FORM.
C      TAPE5       - CARD-IMAGE FILE CONTAINING THE INPUT VALUES AND THE
C                    PARAMETERS FOR A GIVEN RUN OF THIS PROGRAM. MOST OF
C                    THIS FILE IS CONTAINED IN THE RUNVAL NAMELIST.
C      TAPE7       - DATA FILE CONTAINING THE VALUES FOR THE QDUST
C                    ARRAY (DUST TO AIR MIXING RATIO) IN
C                    UNFORMATTED FORM. THE EXACT USAGE OF THIS FILE IS
C                    BE DETERMINED.
C      TAPE9       - DATA FILE CONTAINING THE VALUES FOR THE FLUXNET
C                    ARRAY (SOLAR NET FLUXES) IN UNFORMATTED FORM.
C      TAPE11      - HISTORY TAPE DATA FILE CONTAINING THE ATMOSPHERIC
C                    VARIABLES FOR EACH GRID POINT IN UNFORMATTED FORM.
C                    A HISTORY TAPE FILE IS USED TO INITIALIZE THE SET
C                    OF VALUES FOR EACH GRID POINT IF A COLD START OF
C                    THE PROGRAM IS NOT DONE. ADDITIONAL HISTORY TAPES
C                    MAY BE USED IN INPUT, BUT UNDER NORMAL USAGE, NO
C                    MORE THAN ONE HISTORY TAPE IS USED. (SEE THE
C                    SUBROUTINE NAMED 'INPUT' FOR MORE INFORMATION.)
C  PROGRAM OUTPUT
C      TAPE6       - A PRINT-IMAGE FORMATTED FILE CONTAINING SELECT
C                    DETAILS OF THE EXECUTION OF THE MODEL FOR A RUN.
C                    (FOR DEBUGGING, AND CASUAL ANALYSIS OF RUNS.)
C      TAPE10      - A PRINT-IMAGE FORMATTED FILE CONTAINING MESSAGES
C                    GIVING AN OVERVIEW OF THE EXECUTION OF THE MODEL
C                    FOR A RUN.
C      TAPE11 TO   - HISTORY TAPE FILES. THE PRIMARY (AND MOST
C       TAPE15       VOLUMINOUS) OUTPUT OF THIS PROGRAM IS WRITTEN TO
C                    HISTORY TAPE FILES. FROM ONE TO FIVE HISTORY TAPE
C                    FILES CAN BE GENERATED IN A SINGLE RUN OF THIS
C                    PROGRAM DEPENDING ON THE LENGTH OF THE RUN. EACH
C                    HISTORY TAPE FILE PRODUCED IS STAGED OUT TO A
C                    CORRESPONDING MAGNETIC TAPE. THE HISTORY TAPE
C                    FILES CONTAIN THE ATMOSPHERIC VARIABLES FOR EACH
C                    GRID POINT IN UNFORMATTED FORM.
C
C  SUBROUTINES CALLED:
C      INPUT, OUTPUT, STEP, COMP3, GMP, SDET, HISTORY
C#####################################################################
C                     EXPLANATION OF SELECT VARIABLES
C
C        TAUID...STARTING DAY
C        TAUIH...STARTING HOUR
C        TAUO.....OUTPUT INTERVAL                                      *
C        TAUD.....INTERVAL FOR INCREMENTING SOLAR DECLINATION (24 HRS) *
C        TAUH.....INTERVAL  FOR UPDATING THE HISTORY TAPE (IN HOURS)   *
C        TAUE.....LENGTH OF RUN IN HOURS. THEN CONVERTED TO ENDING TIME*
C
C      CONTROL FLAGS
C
C        KEY1.....TRUE WHEN RUN IS FINISHED.
C        KEY2.....TRUE WHEN A MESSAGE IS DESIRED FOR EACH COMPLETED TIME
C                  STEP.
C        KEY8.....DUMMY FLAG FOR SUBROUTINE OUTPUT.
C        KEY9.....DUMMY FLAG FOR SUBROUTINE OUTPUT.
C        KEY13....TRUE TO SUPRESS HISTORY OUTPUT EXCEPT FOR LAST TIME
C                 STEP.
C        KEY15....TRUE WHEN A MESSAGE IS DESIRED EACH TIME OUTPUT IS
C                 WRITTEN TO THE HISTORY TAPE FILE.
C######################################################################

      use grid_h
      use defines_h
      use constants_h, only: kapa
      use radinc_h
      use radcommon_h
      use standard_h
      use dtcommon_h
      use cldcommon_h
      use comp3cmn_h, only: tauts, taute, volno

      implicit none

      LOGICAL  EVENT, OUT

      integer :: i, j, l
      real*8  :: xtau, presm, testptrop

C     EVENT IS A LOGICAL-VALUED FUNCTION. EVENT(X) IS TRUE ONCE
C     EVERY X HOURS.

c     EVENT(XTAU) = MOD(NSTEP,IFIX(XTAU*DAY/(DT*ROTPER) +.1)).EQ.0
      EVENT(XTAU) = MOD(NSTEP,INT(XTAU*DAY/(DT*ROTPER) +.1)).EQ.0

C     THE FOLLOWING ROW OF = IS MY SYMBOL FOR THE START OF EXECUTABLE
C     CODE IN A PROGRAM MODULE.

C#======================================================================

C     Initialize control flags to their default values

      KEY1  = .FALSE.
      KEY2  = .FALSE.
      KEY8  = .FALSE.
      KEY9  = .FALSE.
      KEY13 = .FALSE.
      KEY15 = .FALSE.

      RESTRT  = .TRUE.
      GOODINP = .TRUE.

C  Set up the radiation code stuff

      CALL RADSETUP

C  Back to the original GCM code
 
      CALL INPUT

      testptrop = ptrop/2.0D0
      if(testptrop.lt.ptop) then
        write(6,'("PTROP must be larger than ",1pe10.3)') 2*ptop
        stop
      end if

C     ABORT PROGRAM IF BAD INPUT.

      IF ( .NOT.  GOODINP  )   GOTO 8000

C     Input is OK, proceed with the program.

      NSTEP  = TAU*DAY/(DT*ROTPER) + .1

      CALL OUTPUT(TAU,ROTPER,RESTRT,IDAY,KEY8,KEY9)

      RESTRT = .FALSE.
C
C     MAIN COMPUTATIONAL CONTROL (TIME STEP) LOOP
C
C     Start the next time step.

200   NSTEP  = NSTEP + 1
 
      TAU    = FLOAT(NSTEP)*ABS(DT)*ROTPER/ DAY   +1.E-3
      TOFDAY = MOD(TAU,ROTPER)   ! Use generic mod
c     TOFDAY = AMOD(TAU,ROTPER)
      KEY8   = KEY8 .OR. EVENT(TAUO)
      KEY9   = KEY9 .OR. EVENT(TAUO2)
      OUT    = (MOD(NSTEP,NC3).EQ.0).AND.( KEY8 .OR. KEY9 )

C     Compute values for the major variables for this time step.

C Turn temperature into theta

      DO I = 1,L_ISIZE
        DO J = 1,L_JSIZE
          DO L =1,L_LAYERS
            PRESM = (SIGMA(2*L+2)*P(J,I))+PTROP
            T(J,I,L) = T(J,I,L)*((PSF/PRESM)**KAPA)
          END DO
        END DO
      END DO

      CALL NEWSTEP
       
C Turn theta into temperature

      DO I = 1,L_ISIZE
        DO J = 1,L_JSIZE
          DO L =1,L_LAYERS
            PRESM = (SIGMA(2*L+2)*P(J,I))+PTROP
            T(J,I,L) = T(J,I,L)*((PRESM/PSF)**KAPA)
          END DO
        END DO
      END DO

      CALL COMP3

C     VARIOUS CHECKING AND HISTORY OPTIONS

      IF(EVENT(24.0D0))  CALL GMP

      IF(EVENT(TAUD)) CALL SDET

      IDAY = TAU/ROTPER

C     IF TIME FOR END OF RUN IS HERE, GO ON UNTIL MOD(NSTEP,NC3).EQ.0
C     (A GOOD STOPPING PLACE) AND WRITE FINAL VALUES TO HISTORY TAPE.

      IF((TAU.GE.TAUE).AND.(MOD(NSTEP,NC3).EQ.0)) KEY1 = .TRUE.
      IF(OUT .OR. KEY1) CALL OUTPUT(TAU,ROTPER,RESTRT,IDAY,KEY8,KEY9)
      IF((EVENT(TAUH).AND.(.NOT. KEY13)).OR.KEY1) CALL HISTORY

      IF(.NOT.KEY1) GO TO 200

C     ############################### END OF MAIN TIME STEP LOOP###
C
C     TAU GREATER THAN TAUE  (TIME TO GO HOME)
C
      WRITE(MTP,9715) IDAY,TOFDAY
      ENDFILE MTP

8000  CONTINUE
      write(6,'(32x,"42")')
      STOP

9200  FORMAT ('1',1X,'MARS GENERAL CIRCULATION MODEL NOW RUNNING')
9715  FORMAT (' ','WMSG036 HAS STOPPED AT DAY ', I4,
     *                   ' / HOUR ',  2 F7.3)
9900  FORMAT (' ',5X, 'STARTING STEP.   NSTEP =', I4,
     *                   '.  ELAPSED EXECUTION TIME =', F7.2 )
9910  FORMAT (' ', 5X, 'STARTING COMP3.  NSTEP =', I4,
     *                   '.  ELAPSED EXECUTION TIME =', F7.2 )

      END
