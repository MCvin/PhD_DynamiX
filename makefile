F90=f90
FFLAGS= -O3 -C
LDFLAGS= 
LIBPATH=/opt/NAGnl/flsol19da/libnag.a
SRC=FPA.f90 photon.f90 telescope.f90 senseurs.f90 outils.f90 lecture.f90
MODPATH=
OBJ=$(SRC:.f90=.o)

all:DynamiX


DynamiX: simu.o
	$(F90) $(FFLAGS) -o $@ simu.o $(OBJ) $(LIBPATH)
	
simu.o: FPA.o photon.o senseurs.o telescope.o
	$(F90) $(FFLAGS) -o $@ -c simu.f90 
	
FPA.o: outils.o lecture.o
	$(F90) $(FFLAGS) -o $@ -c FPA.f90
	
photon.o: outils.o lecture.o
	$(F90) $(FFLAGS) -o $@ -c photon.f90
	
telescope.o: outils.o lecture.o
	$(F90) $(FFLAGS) -o $@ -c telescope.f90
	
senseurs.o: outils.o lecture.o
	$(F90) $(FFLAGS) -o $@ -c senseurs.f90
	
outils.o:
	$(F90) $(FFLAGS) -o $@ -c outils.f90
	
lecture.o:
	$(F90) $(FFLAGS) -o $@ -c lecture.f90
	
clean:
	rm *.o
