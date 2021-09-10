      subroutine getvdust(Ls,nvdust,vdustls,tautotji)

!     Legacy Mars GCM v24
!     Mars Climate Modeling Center
!     NASA Ames Research Center

C
C  Get the opacity map for the requested Ls, return it in TAUTOTJI.
C  In practice, if the Ls values for each map are stored in VDUSTLS.
C  The requested Ls will fall between two of these values, so the
C  returned opacity, TAUTOTJI will be an interpolation between the
C  two opacity maps that bracket the requested Ls.  If the requested
C  Ls is less than VDUSTLS(1), the first opacity map is returned.
C  If Ls > VDUSTLS(NVDUST), the last opacity map is returned.
C
C  It is assumed that VDUSTLS(1) < VDUSTLS(2) < VDUSTLS(3) . . .
C
C  VDUST UPDATES 7/19/01
C
C  INPUT:
C     Ls       - The requested seasonal date
C     NVDUST   - The number of elements in the VDUSTLS array
C     VDUSTLS  - Array of Ls values, which correspond to the Ls
C                of each opacity map in the file (accessed by unit 85)
C
C  OUTPUT:
C     TAUTOTJI - 2D array of the opacity map (Dust opacity at the
C                reference pressure - usually 6.1 mbar)
C
C======================================================================C

      use grid_h
      implicit none

      integer J, I, N, NVDUST, IOS
      real*8 vdustls(1000), TAUTOTJI(L_J,L_I), Ls, RATIO, lsin, lsu
      real*8 opacity1(L_J+1,L_I), opacity2(L_J+1,L_I)
      real*8 opacity(L_J+1,L_I)

!  The opacity maps are now 9-micron - VIS2IR converts to visible
!
      real*8 scale2vis

C======================================================================C

C  Modified 3/6/04 to interpolate through the 360->0 point if we're
C  doing an annual run.  This is assumed to be the case if the Ls of
C  end opacity maps is within 10 degrees of 0/360.
C

      if((Ls.lt.VDUSTLS(1) .or. Ls.gt.VDUSTLS(NVDUST))) then

          read(85,rec=NVDUST) Lsin, opacity1
          read(85,rec=1) Lsin, opacity2

          if(Ls.lt.10.0) then
            lsu = ls + 360.0
          else
            lsu = ls
          end if

          RATIO = (Lsu-VDUSTLS(NVDUST))/
     *            (VDUSTLS(1) + 360.0 -VDUSTLS(NVDUST))

C  C-grid J goes form 1 to JM-1, and the J=1 c-grid PI point
C  is the J=2 b-grid point.

          do J=1,L_J-1
            do I=1,L_I
              TAUTOTJI(J,I) = OPACITY1(J+1,I) + RATIO*
     *                         (OPACITY2(J+1,I) - OPACITY1(J+1,I))
            end do
          end do

      else

C  Interpolate the value

        do n=1,NVDUST-1
          if(Ls.ge.VDUSTLS(n) .and. Ls.lt.VDUSTLS(n+1)) then
            read(85,rec=N) Lsin, opacity1
            read(85,rec=N+1) Lsin, opacity2
            RATIO = (Ls-VDUSTLS(N))/(VDUSTLS(N+1)-VDUSTLS(N))

C  C-grid J goes form 1 to JM-1, and the J=1 c-grid PI point
C  is the J=2 b-grid point.

            do J=1,L_J-1
              do I=1,L_I
                TAUTOTJI(J,I) = OPACITY1(J+1,I) + RATIO*
     *                         (OPACITY2(J+1,I) - OPACITY1(J+1,I))
              end do
            end do

            exit
          end if
        end do

      endif

!  Newer TES opacity maps are 9-micon, older ones are visible
!  Scale the IR maps to visible (2.75 is current TES estimate)
!  Set scale2vis = 2.75 if TES 9-micron maps are used, comment out
!  otherwise.

!  limit dust opacity at 75N, 80N, and 85N to 0.05 (9-micron at 6.1
!  mbar)

!      do I=1,L_I
!        tautotji(35,i) = min(0.025,tautotji(35,I))
!        tautotji(34,i) = min(0.025,tautotji(34,I))
!        tautotji(33,i) = min(0.050,tautotji(33,I))
!      end do


!      scale2vis = 2.75
      scale2vis = 3.67

      do J=1,L_J-1
        do I=1,L_I
          TAUTOTJI(J,I) = scale2vis*TAUTOTJI(J,I)
        end do
      end do

      return
      end
