      subroutine solarza(I,IM,sind,cosd,sinl,cosl,tofday,day,rotper,
     *                   ndt,dt,acosz)

!     Legacy Mars GCM v24
!     Mars Climate Modeling Center
!     NASA Ames Research Center

!     Compute the cosine of the solar zenith angle.  Take into account
!     if the sun rises or sets during the time period under 
!     consideration.

      use constants_h, only: PI

      implicit none

      integer :: i, im
      real*8  :: x3, x6, sind, sinl, cosd, cosl, tofday, day, rotper
      real*8  :: dt, t1, t2, cosz1, cosz2, ndt, x4, x5, acosz, tr1, trr

!======================================================================!

      X5 = 2.0*PI/DAY
      X3 = SIND*SINL
      X6 = 2.0*PI*FLOAT(I-1)/FLOAT(IM)

C     COSINE OF THE SOLAR ZENITH ANGLE ARE DETERMINED ACCORDING TO THE
C     TYPE OF 'PI' POINT.

C     FOR ALL GRID POINTS (WE NO LONGER HAVE POLE POINTS)

      X4    = COSD*COSL
      T1    = (TOFDAY*DAY/ROTPER)-0.5*DT
      T2    = T1 + NDT
      COSZ1 = X3 + X4*COS(X5*T1+X6)
      COSZ2 = X3 + X4*COS(X5*T2+X6)

      IF(COSZ1.GE.0.0.AND.COSZ2.GE.0.0) THEN

C       SUN ABOVE THE HORIZON FOR THE ENTIRE PERIOD.

        ACOSZ = X3+X4*(SIN(X5*T2+X6)-SIN(X5*T1+X6))/(NDT*X5)

      ELSE IF(COSZ1.LE.0.0.AND.COSZ2.LE.0.0) THEN

C       SUN BELOW THE HORIZON FOR THE ENTIRE PERIOD.

        ACOSZ = 0.0

      ELSE

C       SUN RISES OR SETS DURING THE TIME INTERVAL.

        TR1 = ACOS(-X3/X4)
        TRR = (TR1-X6)/X5

C       GET 'TRR' THAT IS BETWEEN T1 AND T2.

        IF(TRR.GE.T1.AND.TRR.LE.T2) GOTO 999
        IF(TRR.GT.T2) GOTO 998
        TRR = TRR+DAY
        IF(TRR.GE.T1.AND.TRR.LE.T2) GOTO 999
        TRR = TRR+DAY
        IF(TRR.GE.T1.AND.TRR.LE.T2) GOTO 999

  998 CONTINUE

        TR1 = 2.0*PI-TR1
        TRR = (TR1-X6)/X5

        IF(TRR.GE.T1.AND.TRR.LE.T2) GOTO 999
        IF(TRR.LT.T1) TRR = TRR+DAY
        IF(TRR.GT.T2) TRR = TRR-DAY

  999 CONTINUE

        IF(COSZ1.LT.0.0.AND.COSZ2.GT.0.0) THEN

C         SUN RISES AFTER TIME T1.

          IF(TRR.LT.T1.OR.TRR.GT.T2) THEN
            TRR = T2
          END IF

          ACOSZ = (X3*(T2-TRR)+X4*(SIN(X5*T2+X6)-
     *                        SIN(X5*TRR+X6))/X5)/NDT

        ELSE

C         SUN SETS BEFORE TIME T2.

          IF(TRR.LT.T1.OR.TRR.GT.T2) THEN
            TRR = T1
          END IF

          ACOSZ=(X3*(TRR-T1)+X4*(SIN(X5*TRR+X6)-SIN(X5*T1+X6))/X5)/NDT

        END IF

      END IF

      return
      end
