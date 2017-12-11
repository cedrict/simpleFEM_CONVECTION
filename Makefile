.SUFFIXES:.out .o .s .c .F .f .f90 .e .r .y .yr .ye .l .p .sh .csh .h

include Makefile.machine

OBJECTS2D = output_for_paraview.o int_to_char.o simplefem.o

.f.o:
	$(F90) $(FLAGS) $(INCLUDE) $*.f
.f90.o:
	$(F90) $(FLAGS) $(INCLUDE) $*.f90

code:	$(OBJECTS2D)
	$(F90) $(OPTIONS) $(OBJECTS2D) $(LIBS) -o simplefem

