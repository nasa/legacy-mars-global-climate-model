!     Legacy Mars GCM v24
!     Mars Climate Modeling Center
!     NASA Ames Research Center

EXE      = .
#F90_COMP = ifort
#F_OPTS   = -c -O2
F90_COMP= gfortran
F_OPTS  = -c -O3 -finit-local-zero -frecord-marker=4

MAKEFILE = Makefile
RUNNAME  = gcm2.3

INCL     = modules.f90

SRC =    main.f input.f output.f comp3.f gmp.f sdet.f history.f\
         emiss.f magfac.f init1.f insdet.f\
         tempgr.f cmp3set.f gridvel.f potemp1.f potemp2.f convect.f\
         cmp3out.f geopotn.f coldair.f getvdust.f jsrchgt.f90\
         jsrchge.f90 mhistv.f mhisth.f gcmlog.f\
         m4read.f bndcond.f descale.f eddycoef.f matrix.f scldef.f\
         solve.f newpbl.f initpbl.f interpdust.f nextls.f dustprofile.f\
         funcd.f newtg.f newstep.f dycore.f filters.f\
         fallvel.f isource.f ssource.f initdt.f initcld.f sizedist.f\
         radsetup.f90 optci.f90 optcv.f90 sfluxi.f90 sfluxv.f90\
         setspi.f90 setspv.f90\
         setrad.f90 tpindex.f90 dsolver.f90 gfluxi.f90 gfluxv.f90\
         dtridgl.f90\
         laginterp.f90 lagrange.f90 filltaucum.f fillpt.f90 dsolflux.f\
         growthrate.f microphys.f opt_cld.f ddevil.f opt_dst.f\
         getdetau.f90\
         ini_optcld.f ini_optdst.f settozero.f solarza.f modules.f90

F90_OBJ = $(SRC:.f=.o)

$(RUNNAME): modules.o $(F90_OBJ) $(INCL) $(MAKEFILE)
	$(F90_COMP) -o $(EXE)/$(RUNNAME) $(F90_OBJ)

modules.o: modules.f90
	$(F90_COMP) $(F_OPTS) $*.f90

.f.o : modules.o
	$(F90_COMP) $(F_OPTS) $*.f

clean:
	rm -f *.o *.mod $(RUNNAME)

$(F90_OBJ): $(MAKEFILE) $($@:.o=.f)
