SRCDIR=src
VPATH=$(SRCDIR)

INC = -I/home/rouckas/include
LIB = -lumfpack -lm -lamd -llapack -lpthread -lboost_filesystem -lboost_system
DEFINES = -DAMD64

CC = g++
CFLAGS =  -Wall
CFLAGS = -g -pg -Wall

OPTIM=-march=core2 -mfpmath=sse -O3 -ffast-math -msse
OPTIM=-march=core2 -mfpmath=sse -O3 -ffast-math -msse
OPENMP=
OPENMP=-fopenmp

.c.o:
	$(CC) $(CFLAGS) $(OPTIM) $(DEFINES) -c $<
.cpp.o:
	$(CC) $(CFLAGS) $(OPTIM) $(DEFINES) $(INC) -c $<


SOURCES = Makefile random.cpp tabulate.cpp input.cpp pic.cpp param.cpp output.cpp matrix.cpp timer.hpp

OBJ = gnuplot_i.o histogram.o particles.o fields.o speclist.o input.o param.o parser.o
OBJ3D = fields3d.o species3d.o

plasma2d: test.cpp $(SOURCES) $(OBJ)
	$(CC) $(CFLAGS) $(DEFINES) -o $@ $(SRCDIR)/test.cpp $(OBJ) $(OPTIM) $(INC) $(LIB) 
plasma3d: plasma3d.cpp $(SOURCES) $(OBJ) $(OBJ3D)
	$(CC) $(CFLAGS) $(DEFINES) -o $@ $(SRCDIR)/plasma3d.cpp $(OBJ) $(OBJ3D) $(OPTIM) $(INC) $(LIB)

penning: penning.cpp $(SOURCES) $(OBJ)
	$(CC) $(CFLAGS) $(DEFINES) -o $@ $(SRCDIR)/penning.cpp $(OBJ) $(OPTIM) $(INC) $(LIB) 
CRDS: CRDS.cpp $(SOURCES) $(OBJ)
	$(CC) $(CFLAGS) $(DEFINES) -o $@ $(SRCDIR)/CRDS.cpp $(OBJ) $(OPTIM) $(INC) $(LIB)
CRDS1D: CRDS1D.cpp $(SOURCES) $(OBJ)
	$(CC) $(CFLAGS) $(DEFINES) -o $@ $(SRCDIR)/CRDS1D.cpp $(OBJ) $(OPTIM) $(INC) $(LIB)
MAC: MAC.cpp $(SOURCES) $(OBJ)
	$(CC) $(CFLAGS) $(DEFINES) -o $@ $(SRCDIR)/MAC.cpp $(OBJ) $(OPTIM) $(INC) $(LIB)
#test_scatter_speed: test_scatter_speed.cpp particles.cpp gnuplot_i.c argon.cpp elon.cpp Makefile random.cpp tabulate.cpp input.cpp elonO2.cpp
#	$(CC) $(CFLAGS) -o $@ test_scatter_speed.cpp gnuplot_i.c histogram.cpp $(OPTIM) $(INC) $(LIB) 
#test_flux: test_flux.cpp random.c random.h
#	$(CC) $(CFLAGS) $(OPTIM) -o ../$@ test_flux.cpp random.c histogram.cpp
#rot_test: rot_test.cpp random.c random.h
#	$(CC) $(CFLAGS) $(OPTIM) -o ../$@ rot_test.cpp random.c
#maxwell_test: maxwell_test.cpp random.c random.h
#	$(CC) $(CFLAGS) $(OPTIM) -o ../$@ maxwell_test.cpp random.c gnuplot_i.c histogram.cpp
test_random_omp: tests/test_random_omp.cpp src/random.cpp Makefile
	$(CC) $(CFLAGS) $(DEFINES) -o $@ tests/test_random_omp.cpp $(OPTIM) $(OPENMP)
test_particles_omp: tests/test_particles_omp.cpp src/random.cpp Makefile
	$(CC) $(CFLAGS) $(DEFINES) -o $@ tests/test_particles_omp.cpp $(OPTIM) $(OPENMP)
gnuplot_i.o: gnuplot_i.h
histogram.o: histogram.hpp
particles.o: particles.hpp random.cpp mymath.cpp tabulate.cpp util.cpp param.hpp fields.hpp parser.hpp
parser.o: parser.hpp
fields.o: fields.hpp
speclist.o: speclist.hpp
input.o: input.hpp
param.o: param.hpp
fields3d.o: fields3d.hpp
species3d.o: species3d.hpp
clean:
	rm -f *.o
