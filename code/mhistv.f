      SUBROUTINE MHISTV(TAU,P,T,U,V,GT,CO2ICE,TAUSURF,STRESSX,STRESSY,
     *                  DXYP,
     *                  TSTRAT,SSUN,QTRACE,QCOND,STEMP,NC3,NCYCLE,PTROP,
     *                  vpout,rsdist,tofday,psf,tautot,prtau,sind,
     *                  fuptopv,fdntopv,fupsurfv,fdnsurfv,fuptopir,
     *                  fupsurfir,fdnsurfir,srfupflx,srfdnflx,
     *                  tauref3d,dheat,geot)

!     Legacy Mars GCM v24
!     Mars Climate Modeling Center
!     NASA Ames Research Center

C     Mintz HISTory Variables - WRITE out the history variables
C     for each history time step in real*4 format.  A special
C     subroutine will WRITE out the Mintz HISTory Header record.
C
C     Convert GCM history tapes from ieee real*8 to ieee real*4
C     for output.  A special startup file with real*8 will also be
C     written in HISTORY every time a new history "tape" is started.
C
C     Original starting code was sgi2sun.f
C
C     Version 0.1
C     Started Feb. 3, 1994
C

      use grid_h
      use defines_h
      use dtcommon_h
      use cldcommon_h
      use standard_h, only: npcflag, surfalb

      implicit none

      integer, parameter :: IBLKMX = 161

      real*8  :: area

C----------------------------------------------------------------------#

      REAL*8 TAU, PTROP

      REAL*8 TAUSURF(L_J,L_I), QTRACE(L_J,L_I,L_LAYERS,NTRACE)
      REAL*8 DHEAT(L_J,L_I,L_LAYERS)
      REAL*8 GEOT(L_J,L_I,L_LAYERS)
      REAL*8 QCOND(L_J,L_I,NTRACE)
      real*8 stemp(L_J,L_I,2*NL+1)
      real*8 vpout,rsdist,tofday,psf,tautot,prtau,sind,dxyp(L_J)

      REAL*4 TAUSURF4(L_J,L_I), QTRACE4(L_J,L_I,L_LAYERS,NTRACE)
      REAL*4 DHEAT4(L_J,L_I,L_LAYERS)
      REAL*4 GEOT4(L_J,L_I,L_LAYERS)
      REAL*4 QCOND4(L_J,L_I,NTRACE)
      real*4 stemp4(L_J,L_I,NL)

      REAL*8 P(L_J,L_I), GT(L_J,L_I), CO2ICE(L_J,L_I)
      REAL*8 T(L_J,L_I,L_LAYERS)
      REAL*8 U(L_J,L_I,L_LAYERS),V(L_J,L_I,L_LAYERS)
      REAL*8 TSTRAT(L_J,L_I), SSUN(L_J,L_I)
      REAL*8 STRESSX(L_J,L_I), STRESSY(L_J,L_I)

      REAL*4 P4(L_J,L_I),GT4(L_J,L_I)
      REAL*4 T4(L_J,L_I,L_LAYERS)
      REAL*4 U4(L_J,L_I,L_LAYERS),V4(L_J,L_I,L_LAYERS)
      REAL*4 CO2ICE4(L_J,L_I)
      REAL*4 STRESSX4(L_J,L_I), STRESSY4(L_J,L_I)
      REAL*4 TSTRAT4(L_J,L_I), SSUN4(L_J,L_I)

      REAL*8 RUNNUM
      REAL*4 RUNNUM4
      REAL*4 RPTAU4
      real*4 gasp4

      REAL*8 SRFUPFLX(L_J,L_I,NTRACE)
      REAL*8 SRFDNFLX(L_J,L_I,NTRACE)
      REAL*8 TAUREF3D(L_J,L_I,L_LEVELS)
      REAL*4 SRFUPFLX4(L_J,L_I,NTRACE)
      REAL*4 SRFDNFLX4(L_J,L_I,NTRACE)
      REAL*4 TAUREF3D4(L_J,L_I,L_LEVELS)

C  Dust tracer variables

      REAL*4 DPDEN4, DRAD4(NDP), DR4(NDP), QEXT4(NDP), GTAU4
      REAL*4 SSAREA4, FLUX4(NDP), DSD4, THRESHOLD4

      INTEGER IBLKCT

!  Implicit none

      integer j, i, l, n, nc3, ncycle, k
      real*4  tau4, vpout4, rsdist4, tofday4, psf4, ptrop4, tautot4
      real*4  sind4

!  Radiation variables

      real*4 fuptopv(L_J,L_I), fdntopv(L_J,L_I), fupsurfv(L_J,L_I),
     *       fdnsurfv(L_J,L_I), fuptopir(L_J,L_I), fupsurfir(L_J,L_I),
     *       fdnsurfir(L_J,L_I)

C======================================================================C

      DO J=1,L_JSIZE-1
        DO I=1,L_ISIZE
          TSTRAT4(J,I)   = TSTRAT(J,I)
          TAUSURF4(J,I)  = TAUSURF(J,I)
          SSUN4(J,I)     = SSUN(J,I)
          STRESSX4(J,I)  = STRESSX(J,I)
          STRESSY4(J,I)  = STRESSY(J,I)
          GT4(J,I)       = GT(J,I)
          CO2ICE4(J,I)   = CO2ICE(J,I)
          DO N=1,NTRACE
                SRFUPFLX4(J,I,N) = SRFUPFLX(J,I,N)
                SRFDNFLX4(J,I,N) = SRFDNFLX(J,I,N)
          ENDDO
          DO K=1,L_LEVELS
                TAUREF3D4(J,I,K)=TAUREF3D(J,I,K)
          ENDDO
          DO L=1,NL
            STEMP4(J,I,L) = stemp(J,I,2*L)
          end do
        END DO
      END DO

      DO J=1,L_JSIZE-1
        DO I=1,L_ISIZE
          DO L=1,L_LAYERS
            T4(J,I,L) = T(J,I,L)
            U4(J,I,L) = U(J,I,L)
            V4(J,I,L) = V(J,I,L)
            DHEAT4(J,I,L) = DHEAT(J,I,L)
            GEOT4(J,I,L) = GEOT(J,I,L)
            DO N=1,NTRACE
              if(QTRACE(J,I,L,N).gt.0.0 .and. 
     *                           QTRACE(J,I,L,N).lt.1.0E-30) then
                QTRACE4(J,I,L,N) = 1.0E-30
              else
                QTRACE4(J,I,L,N) = QTRACE(J,I,L,N)
              end if
            END DO
          END DO
          P4(J,I) = P(J,I)

          DO N=1,NTRACE
            if(QCOND(J,I,N).gt.0.0 .and. QCOND(J,I,N).lt.1.0E-30) then
              QCOND4(J,I,N) = 1.0E-30
            elseif
     *       (QCOND(J,I,N).lt.0.0 .and. QCOND(J,I,N).gt.-1.0E-30) then
              QCOND4(J,I,N) = -1.0E-30
            else
              QCOND4(J,I,N) = QCOND(J,I,N)
            end if
          END DO
          if((.not.npcflag(j,i)) .and. qcond4(J,I,iMa_vap).lt.0.0) then
            qcond4(J,I,iMa_vap) = 1.0E-30
          end if
        ENDDO
      ENDDO

      gasp4 = 0.0
      area = 0.0
      do j=1,L_J-1
        do I=1,L_I
          area  = area + dxyp(j)
          gasp4 = gasp4 + (P(J,I)+ptrop)*dxyp(J)
        end do
      end do

      gasp4 = gasp4/area

C  Fill the real*4 dust tracer variables

      TAU4    = tau
      VPOUT4  = vpout
      RSDIST4 = rsdist
      TOFDAY4 = tofday
      PSF4    = psf
      PTROP4  = ptrop
      TAUTOT4 = tautot
      RPTAU4  = prtau
      SIND4   = sind
      
      WRITE(11) TAU4, VPOUT4, RSDIST4, TOFDAY4, PSF4, PTROP4, TAUTOT4,
     *          RPTAU4, SIND4, GASP4
      WRITE(11) NC3, NCYCLE

      WRITE(11) P4
      WRITE(11) T4
      WRITE(11) U4
      WRITE(11) V4
      WRITE(11) GT4
      WRITE(11) CO2ICE4
      WRITE(11) STRESSX4
      WRITE(11) STRESSY4
      WRITE(11) TSTRAT4
      WRITE(11) TAUSURF4
      WRITE(11) SSUN4
      WRITE(11) QTRACE4
      WRITE(11) QCOND4
      write(11) STEMP4
      write(11) fuptopv, fdntopv, fupsurfv, fdnsurfv
      write(11) fuptopir, fupsurfir, fdnsurfir
      write(11) surfalb
      write(11) dheat4
      write(11) geot4

      print*,srfupflx4(12,12,1),srfdnflx4(12,12,1)
      print*,srfupflx4(35,1,6),srfdnflx4(35,1,6)
      print*,'tauref3d4(17,1,35)',tauref3d4(17,1,35)

      write(45) tau4, vpout4, tofday4, srfupflx4, srfdnflx4,
     *          tauref3d4

 
      RETURN
      END
