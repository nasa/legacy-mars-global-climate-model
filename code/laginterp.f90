      subroutine laginterp(pgref,pint,co2i,co2v,fzeroi,fzerov)

!     Legacy Mars GCM v24
!     Mars Climate Modeling Center
!     NASA Ames Research Center

!  Lagrange interpolation (linear in log pressure) of the CO2 
!  k-coefficients in the pressure domain.  Subsequent use of these
!  values will use a simple linear interpolation in pressure.

      use grid_h
      use radinc_h

      implicit none

      integer :: n, nt, np, nh, ng, nw, m, i
      real*8  :: co2i8(L_NTREF,L_NPREF,L_REFH2O,L_NSPECTI,L_NGAUSS)
      real*8  :: co2v8(L_NTREF,L_NPREF,L_REFH2O,L_NSPECTV,L_NGAUSS)
      real*8  :: pgref(L_NPREF)

      real*8  :: co2i(L_NTREF,L_PINT,L_REFH2O,L_NSPECTI,L_NGAUSS)
      real*8  :: co2v(L_NTREF,L_PINT,L_REFH2O,L_NSPECTV,L_NGAUSS)

      real*8  :: fzeroi(L_NSPECTI)
      real*8  :: fzerov(L_NSPECTV)
 
      real*8  :: x, xi(4), yi(4), ans
      real*8  :: pint(L_PINT), pref(L_NPREF), p

      real*8  :: pin(L_PINT) =  [                                      &
                           -6.0D0, -5.8D0, -5.6D0, -5.4D0, -5.2D0,     &
                           -5.0D0, -4.8D0, -4.6D0, -4.4D0, -4.2D0,     &
                           -4.0D0, -3.8D0, -3.6D0, -3.4D0, -3.2D0,     &
                           -3.0D0, -2.8D0, -2.6D0, -2.4D0, -2.2D0,     &
                           -2.0D0, -1.8D0, -1.6D0, -1.4D0, -1.2D0,     &
                           -1.0D0, -0.8D0, -0.6D0, -0.4D0, -0.2D0,     &
                            0.0D0,  0.2D0,  0.4D0,  0.6D0,  0.8D0,     &
                            1.0D0,  1.2D0,  1.4D0,  1.6D0,  1.8D0,     &
                            2.0D0,  2.2D0,  2.4D0,  2.6D0,  2.8D0,     &
                            3.0D0,  3.2D0,  3.4D0,  3.6D0,  3.8D0,     &
                            4.0D0                                  ]

!======================================================================!

!  Fill pint for output from this subroutine

      do n=1,L_PINT
        PINT(n) = PIN(n)
      end do

!  Take log of the reference pressures

      do n=1,L_NPREF
        pref(n) = LOG10(PGREF(n))
      end do

!     Get CO2 k coefficients

      open(20,file='data/CO2H2O_V_12_95_INTEL',   &
              form='unformatted')
      read(20) co2v8
      read(20) fzerov
      close(20)

      open(20,file='data/CO2H2O_IR_12_95_INTEL',  &
              form='unformatted')
      read(20) co2i8
      read(20) fzeroi
      close(20)

!  Take Log10 of the values - we interpolate the log10 of the values,
!  not the values themselves.   Smallest value is 1.0E-200.

      do nt=1,L_NTREF
        do np=1,L_NPREF
          do nh=1,L_REFH2O
            do ng = 1,L_NGAUSS

              do nw=1,L_NSPECTV
                if(co2v8(nt,np,nh,nw,ng).gt.1.0d-200) then
                  co2v8(nt,np,nh,nw,ng) = log10(co2v8(nt,np,nh,nw,ng))
                else
                  co2v8(nt,np,nh,nw,ng) = -200.0
                end if
              end do
  
              do nw=1,L_NSPECTI
                if(co2i8(nt,np,nh,nw,ng).gt.1.0d-200) then
                  co2i8(nt,np,nh,nw,ng) = log10(co2i8(nt,np,nh,nw,ng))
                else
                  co2i8(nt,np,nh,nw,ng) = -200.0
                end if
              end do
      
            end do
          end do
        end do
      end do

!  Interpolate the values:  first the IR

      do nt=1,L_NTREF
        do nh=1,L_REFH2O
        do nw=1,L_NSPECTI
          do ng=1,L_NGAUSS

!  First, the initial interval (P=1e-6 to 1e-5)

            n = 1 
            do m=1,5
              x     = pint(m)
              xi(1) = pref(n)
              xi(2) = pref(n+1)
              xi(3) = pref(n+2)
              xi(4) = pref(n+3)
              yi(1) = co2i8(nt,n,nh,nw,ng)
              yi(2) = co2i8(nt,n+1,nh,nw,ng)
              yi(3) = co2i8(nt,n+2,nh,nw,ng)
              yi(4) = co2i8(nt,n+3,nh,nw,ng)
              call lagrange(x,xi,yi,ans)
              co2i(nt,m,nh,nw,ng) = 10.0**ans
            end do 
 
            do n=2,L_NPREF-2
              do m=1,5
                i     = (n-1)*5+m
                x     = pint(i)
                xi(1) = pref(n-1)
                xi(2) = pref(n)
                xi(3) = pref(n+1)
                xi(4) = pref(n+2)
                yi(1) = co2i8(nt,n-1,nh,nw,ng)
                yi(2) = co2i8(nt,n,nh,nw,ng)
                yi(3) = co2i8(nt,n+1,nh,nw,ng)
                yi(4) = co2i8(nt,n+2,nh,nw,ng)
                call lagrange(x,xi,yi,ans)
                co2i(nt,i,nh,nw,ng) = 10.0**ans
              end do 
            end do

!  Now, get the last interval (P=1e+3 to 1e+4)

            n = L_NPREF-1
      
            do m=1,5
              i     = (n-1)*5+m
              x     = pint(i)
              xi(1) = pref(n-2)
              xi(2) = pref(n-1)
              xi(3) = pref(n)
              xi(4) = pref(n+1)
              yi(1) = co2i8(nt,n-2,nh,nw,ng)
              yi(2) = co2i8(nt,n-1,nh,nw,ng)
              yi(3) = co2i8(nt,n,nh,nw,ng)
              yi(4) = co2i8(nt,n+1,nh,nw,ng)
              call lagrange(x,xi,yi,ans)
              co2i(nt,i,nh,nw,ng) = 10.0**ans
            end do  

!  Fill the last pressure point

            co2i(nt,L_PINT,nh,nw,ng) = 10.0**co2i8(nt,L_NPREF,nh,nw,ng)

          end do
        end do
        end do
      end do

!  Interpolate the values:  now the visible

      do nt=1,L_NTREF
        do nh=1,L_REFH2O
        do nw=1,L_NSPECTV
          do ng=1,L_NGAUSS

!  First, the initial interval (P=1e-6 to 1e-5)

            n = 1 
            do m=1,5
              x     = pint(m)
              xi(1) = pref(n)
              xi(2) = pref(n+1)
              xi(3) = pref(n+2)
              xi(4) = pref(n+3)
              yi(1) = co2v8(nt,n,nh,nw,ng)
              yi(2) = co2v8(nt,n+1,nh,nw,ng)
              yi(3) = co2v8(nt,n+2,nh,nw,ng)
              yi(4) = co2v8(nt,n+3,nh,nw,ng)
              call lagrange(x,xi,yi,ans)
              co2v(nt,m,nh,nw,ng) = 10.0**ans
            end do 
 
            do n=2,L_NPREF-2
              do m=1,5
                i     = (n-1)*5+m
                x     = pint(i)
                xi(1) = pref(n-1)
                xi(2) = pref(n)
                xi(3) = pref(n+1)
                xi(4) = pref(n+2)
                yi(1) = co2v8(nt,n-1,nh,nw,ng)
                yi(2) = co2v8(nt,n,nh,nw,ng)
                yi(3) = co2v8(nt,n+1,nh,nw,ng)
                yi(4) = co2v8(nt,n+2,nh,nw,ng)
                call lagrange(x,xi,yi,ans)
                co2v(nt,i,nh,nw,ng) = 10.0**ans
              end do 
            end do

!  Now, get the last interval (P=1e+3 to 1e+4)

            n = L_NPREF-1
      
            do m=1,5
              i     = (n-1)*5+m
              x     = pint(i)
              xi(1) = pref(n-2)
              xi(2) = pref(n-1)
              xi(3) = pref(n)
              xi(4) = pref(n+1)
              yi(1) = co2v8(nt,n-2,nh,nw,ng)
              yi(2) = co2v8(nt,n-1,nh,nw,ng)
              yi(3) = co2v8(nt,n,nh,nw,ng)
              yi(4) = co2v8(nt,n+1,nh,nw,ng)
              call lagrange(x,xi,yi,ans)
              co2v(nt,i,nh,nw,ng) = 10.0**ans
            end do  

!  Fill the last pressure point

            co2v(nt,L_PINT,nh,nw,ng) = 10.0**co2v8(nt,L_NPREF,nh,nw,ng)
            
          end do
        end do
        end do
      end do

      return
      end subroutine laginterp
