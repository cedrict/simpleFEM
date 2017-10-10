.SUFFIXES:.out .o .s .c .F .f .f90 .e .r .y .yr .ye .l .p .sh .csh .h

include Makefile.mac

OBJECTS2D = linpack_d.o simplefem.o

.f.o:
	$(F90) $(FLAGS) $(INCLUDE) $*.f
.f90.o:
	$(F90) $(FLAGS) $(INCLUDE) $*.f90

code:	$(OBJECTS2D)
	$(F90) $(OPTIONS) $(OBJECTS2D) $(LIBS) -o simplefem

clean:
	-@rm -f simplefem
	-@rm -f *.o

cleandata:
	-@rm -f OUT/solution_u.dat
	-@rm -f OUT/solution_v.dat

