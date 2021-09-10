      SUBROUTINE  MAGFAC

!     Legacy Mars GCM v24
!     Mars Climate Modeling Center
!     NASA Ames Research Center

C
C  PURPOSE:
C      MAGFAC CALCULATES THE DISTANCE BETWEEN ADJACENT GRID POINTS AND
C      THE AREAS OF THE GRID SQUARES, AND CALCULATES THE CORIOLIS
C      PARAMETERS.
C  MODIFIED FROM ORIGINAL BY
C      STEVE POHORSKY    INFORMATICS     TASK 605    JUL 82
C  FOR
C      JIM POLLACK
C  ENVIRONMENT
C      Cray-2            UNICOS 3.0      FORTRAN
C  REVISION HISTORY
C
C  INPUT PARAMETERS
C      RAD     - RADIUS OF MARS IN METERS.
C      DLAT    - THE LATITUDINAL SPACING BETWEEN GRID POINTS IN RADIANS
C                ( PI/36 ).
C      DLON    - THE LONGITUDINAL SPACING BETWEEN GRID POINTS IN RADIANS
C                ( 2*PI/60 ).
C  OUTPUT PARAMETERS
C      LAT(J) ARRAY   - THE LATITUDE OF THE 'PI' GRID POINTS WITH INDEX
C                       J (IN RADIANS.) 
C      DXV(J) ARRAY   - EAST-WEST DISTANCE IN METERS BETWEEN ADJACENT
C                       'PI' GRID POINTS.
C      DXU(J) ARRAY   - EAST-WEST DISTANCE IN METERS BETWEEN ADJACENT
C                       'U,PI' GRID POINTS.
C      DXYP(J) ARRAY  - THE AREA OF THE GRID SQUARE AROUND A 'PI' GRID
C                       POINT AT LATITUDE INDEX J. (THE GRID SQUARE IS
C                       THE TRAPEZOID FORMED BY THE FOUR NEIGHBORING
C                       'U,V' POINTS. )
C  CALLED BY
C      INPUT
C
      use grid_h
      use defines_h
      use constants_h, only: PI
      use standard_h

      implicit none

C######################################################################
      integer :: j
C#=====================================================================

C The C grid variables are offset in a strange way such that the C PI grid
C is offset from the pole by DY, NOT DY/2, even though PI points are 
C not stored at the pole in the C grid
C
      DO J = 1,L_JSIZE-1
        LAT(J) = (REAL(J)*PI/REAL(L_JSIZE))-(PI/2.0)
      END DO

      DO 420 J=1,L_JSIZE-1
        DXU(J) = RAD*COS(LAT(J))*DLON
  420 CONTINUE

      DO 430 J=2,L_JSIZE-1
        DXV(J) = 0.5*(DXU(J)+DXU(J-1))
  430 CONTINUE

      DXV(1)       = 0.0
      DXV(L_JSIZE) = 0.0

      DO 440 J=2,L_JSIZE-2
        DYV(J) = RAD*DLAT
  440 CONTINUE

      DYV(1)         = RAD*DLAT*1.0
      DYV(L_JSIZE-1) = RAD*DLAT*1.0

      DO 445 J=2,L_JSIZE-2
        DXYP(J) = 0.5*(DXV(J)+DXV(J+1))*DYV(J)
  445 CONTINUE
      DXYP(1) = DXU(1)*DYV(1)
      DXYP(L_JSIZE-1) = DXU(L_JSIZE-1)*DYV(L_JSIZE-1)

      RETURN
      END
