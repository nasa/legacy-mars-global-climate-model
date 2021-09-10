      subroutine initdt

!     Legacy Mars GCM v24
!     Mars Climate Modeling Center
!     NASA Ames Research Center

C
C  Initialize the variables used in the dust tracer scheme:  dustend.f

      use grid_h
      use defines_h
      use standard_h
      use dtcommon_h

      implicit none

      integer :: i,j,n

C======================================================================C

C  Some general constants - fill variables that will be passed to the 
C  dust tracer scheme via dtcommon.h.

      nlay_dt = L_LAYERS
      nlev_dt = 2*L_LAYERS+3
      JM_dt   = L_J
      IM_dt   = L_I
      
      scale_dt = 3.0*(1.38066E-23/1.66054E-27)/44.0

      do J=1,L_J-1
        dxyp_dt(J) = dxyp(J)
      end do

C  Dust particle properties
C  Dust particle density

      dpden_dt = 2.5E+3

C  NDP_DT is the number of dust particle bins.  The associated dust particle
C  radius is drad - dimensioned drad(ndp).  DR_DT(NDP) is the radius bin
C  width, Qext_dt(NDP) is the extinction coefficient for each particle
C  bin.

      ndp_dt     = NDP

C  For now, all dust bins have the same particle size.  They differ only 
C  in their threshold values for lifting in subroutine DUSTEND.

      do n=1,ndp    
       drad_dt(n) = 2.5E-6
       dr_dt(n)   = 5.0E-6
       Qext_dt(n) = 2.938   !QextREF(2) - reference wavelength extinction
      end do

C  NOT_DT is the number of Other Tracer particles carried

      not_dt = L_NOT

C  Calculate surface area of the planet

      area_dt = 0.0
      do J=1,L_J-1
        area_dt = area_dt + dxyp(j)*float(L_I)      
      end do

C  Specified source?  specified = true for specified source, interactive
C  source otherwise

c     specified_dt = .true.
      specified_dt = .false.

C  If a specified source is wanted, set up for that, otherwise set up
C  for interactive dust lifting.

      if(specified_dt .eqv. .TRUE.) then

C  Location of specified source

        J1_dt =  6
        J2_dt =  9
        I1_dt =  1
        I2_dt = 40

C  Calculate the area of the source region

        ssarea_dt = 0.0
        DO I = I1_DT,I2_DT
          DO J=J1_DT,J2_DT
            ssarea_dt = ssarea_dt + dxyp(J)
          END DO
        END DO

        write(6,'("Dust source area = ",1pe10.3)') ssarea_dt

C  Global dust opacity

        gtau_dt = 5.00 

C  Duration, seconds, of source lifting

c       dsd_dt = 88775.0    ! 1 sol
        dsd_dt = 10.0*88775.0  

C  Calculate the mass flux from the surface for the specified region
C  assuming NDP particle sizes, and a given size distribution

        call sizedist(GTAU_DT,DSD_DT,NDP_DT,DRAD_DT,DR_DT,
     *                QEXT_DT,AREA_DT,SSAREA_DT,DPDEN_DT,
     *                FLUX_DT)

      else

C  Interactive source threshold value

        threshold_dt = 2.25E-2

      end if

      return
      end
