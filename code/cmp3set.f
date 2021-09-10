      SUBROUTINE CMP3SET(IM,NC3,DT,IMOVER2,IMOVER4,NDT,UT,VT,TT,
     *                   QTDELTA,FRY, HSOLCO2, HSOLDST, HTH15, HTHOUT, 
     *                   HCONADJ, HTURBO, HRINUM, FCONADJ, FTURBO, 
     *                   FRINUM, FRAYFR, HRAYFR, DISRAY, CO2LAT,QCDEL)

!     Legacy Mars GCM v24
!     Mars Climate Modeling Center
!     NASA Ames Research Center

C  PURPOSE:
C      CMP3SET PERFORMS LOOP INITILIZATION FOR THE LARGE ATMOSPHERIC
C      TEMPERATURE CALCULATION LOOP IN COMP3.
C  AUTHOR
C      STEVE POHORSKY    INFORMATICS     TASK 605    OCT 81
C  FOR
C      JIM POLLACK
C  ENVIRONMENT
C      Cray-2            UNICOS 3.0      FORTRAN
C  REVISION HISTORY
C
C  OUTPUT PARAMETERS
C      UT, VT, & TT ARRAYS, AND VARIABLES IN LOOPVAL COMMON BLOCK
C      ARE INITIALIZED.
C
C  CALLED BY
C      COMP3
C
      use grid_h
      use defines_h

      implicit none
 
C######################################################################
C
      REAL*8 :: NDT

      REAL*8 :: UT(L_JSIZE,L_ISIZE,L_LAYERS)
      REAL*8 :: VT(L_JSIZE,L_ISIZE,L_LAYERS)
      REAL*8 :: TT(L_JSIZE,L_ISIZE,L_LAYERS)
      REAL*8 :: QTDELTA(L_JSIZE,L_ISIZE,L_LAYERS,NTRACE)
      REAL*8 :: QCDEL(L_JSIZE,L_ISIZE,NTRACE)
      real*8 :: FRY(L_JSIZE,L_LAYERS)
      real*8 :: HSOLCO2(L_JSIZE,L_LAYERS)
      real*8 :: HSOLDST(L_JSIZE,L_LAYERS)
      real*8 :: HTH15(L_JSIZE,L_LAYERS)
      real*8 :: HTHOUT(L_JSIZE,L_LAYERS)
      real*8 :: HCONADJ(L_JSIZE,L_LAYERS)
      real*8 :: HTURBO(L_JSIZE,L_LAYERS)
      real*8 :: HRINUM(L_JSIZE,L_LAYERS)
      real*8 :: FCONADJ(L_JSIZE,L_LAYERS)
      real*8 :: FTURBO(L_JSIZE,L_LAYERS)
      real*8 :: FRINUM(L_JSIZE,L_LAYERS)
      real*8 :: FRAYFR(L_JSIZE,L_LAYERS)
      real*8 :: HRAYFR(L_JSIZE,L_LAYERS)
      real*8 :: DISRAY(L_JSIZE,L_LAYERS)
      real*8 :: CO2LAT(L_JSIZE,L_LAYERS)

!     implicit none

      integer :: I, J, L, M, IM, JM, IMOVER2, IMOVER4, NC3
      real*8  :: DT

C=======================================================================

      IMOVER2 = IM/2
      IMOVER4 = IM/4

C     Atmospheric temperature calculation time step

      NDT   = FLOAT(NC3)*DT

C     ZERO LARGE ARRAYS.

      DO L=1,L_LAYERS
        DO J=1,L_JSIZE
          DO I=1,L_ISIZE
            UT(J,I,L) = 0.0
            VT(J,I,L) = 0.0
            TT(J,I,L) = 0.0
            DO M = 1,NTRACE
              QTDELTA(J,I,L,M) = 0.0
            END DO
          END DO
        END DO
      END DO

      DO J=1,L_JSIZE
        DO I=1,L_ISIZE
          DO M = 1,NTRACE
            QCDEL(J,I,M) = 0.0D0
          END DO
        END DO
      END DO

C
C     Set the (J,L) arrays to zero each time through COMP3.
C
      DO 1030 L=1,L_LAYERS
           DO 1010 J=1,L_JSIZE
             FRY(J,L)     = 0.0
             HSOLCO2(J,L) = 0.0
             HSOLDST(J,L) = 0.0
             HTH15(J,L)   = 0.0
             HTHOUT(J,L)  = 0.0
             HCONADJ(J,L) = 0.0
             HTURBO(J,L)  = 0.0
             HRINUM(J,L)  = 0.0
             FCONADJ(J,L) = 0.0
             FTURBO(J,L)  = 0.0
             FRINUM(J,L)  = 0.0
             FRAYFR(J,L)  = 0.0
             HRAYFR(J,L)  = 0.0
             DISRAY(J,L)  = 0.0
             CO2LAT(J,L)  = 0.0
 1010      CONTINUE
 1030 CONTINUE

 3900 CONTINUE
      RETURN
      END
