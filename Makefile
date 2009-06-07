#LIB = -lumfpack -lm -lamd -lgoto_core2p-r1.19 -lgfortran -lboost_filesystem
#LIB = -lumfpack -lm -lamd -lgoto -lpthread -L/opt/UMFPACK/Lib -L/opt/AMD/Lib -L/opt/GotoBLAS
LIB = -lumfpack -lm -lamd  -lpthread -L/opt/UMFPACK/Lib -L/opt/AMD/Lib -L/opt/GotoBLAS
#LIB = -lumfpack -lm -lamd -lgoto_athlon-r1.11

#kevf23 flags
LIB = -lumfpack -lm -lamd -lgoto_northwoodp-r1.19 -lboost_filesystem
INC = -I/usr/local/include/UMFPACK/ -I/usr/local/include/UFconfig/ -I/usr/local/include/AMD/ -I/opt/UMFPACK/Include -I/opt/UFconfig -I/opt/AMD/Include
DEFINES = 

#amd64 flags
INC = -I/opt/UMFPACK/Include -I/opt/UFconfig -I/opt/AMD/Include -I/home/stepan/include
LIB = -lumfpack -lm -lamd -lgoto -lpthread -lboost_filesystem -L/opt/UMFPACK/Lib -L/opt/AMD/Lib -L/opt/GotoBLAS
DEFINES = -DAMD64

#kevf132 flags
INC = -I/opt/UMFPACK/Include -I/opt/UFconfig -I/opt/AMD/Include -I/home/stepan/include
LIB = -lumfpack -lm -lamd -llapack -lpthread -lboost_filesystem -L/opt/UMFPACK/Lib -L/opt/AMD/Lib -L/opt/GotoBLAS
DEFINES = -DAMD64

CC = g++
CFLAGS =  -Wall
CFLAGS = -g -pg -Wall

MODUL=poisson.c random.c
OBJ=gnuplot_i.o
OPTIM=-march=nocona	-mfpmath=sse -O3 -ffast-math -msse
OPTIM=-march=core2	-mfpmath=sse -O3 -ffast-math -msse
#OPTIM=-march=athlon-xp -mfpmath=sse -O3 -ffast-math -msse

.c.o:
	$(CC) $(CFLAGS) -c $<
.cpp.o:
	$(CC) $(CFLAGS) -c $<


SOURCES = particles.cpp argon.cpp elon.cpp Makefile random.cpp tabulate.cpp input.cpp elonO2.cpp pic.cpp oxygen.cpp param.cpp output.cpp fields.cpp argonO2.cpp histogram.cpp matrix.cpp speclist.hpp hydrogen.cpp timer.hpp

OBJ = gnuplot_i.o histogram.o

plasma2d: test.cpp $(SOURCES) $(OBJ)
	$(CC) $(CFLAGS) $(DEFINES) -o $@ test.cpp $(OBJ) $(OPTIM) $(INC) $(LIB) 
penning: penning.cpp $(SOURCES) $(OBJ)
	$(CC) $(CFLAGS) $(DEFINES) -o $@ penning.cpp $(OBJ) $(OPTIM) $(INC) $(LIB) 
MAC: MAC.cpp $(SOURCES) $(OBJ)
	$(CC) $(CFLAGS) $(DEFINES) -o $@ MAC.cpp $(OBJ) $(OPTIM) $(INC) $(LIB)
test_scatter_speed: test_scatter_speed.cpp particles.cpp gnuplot_i.c argon.cpp elon.cpp Makefile random.cpp tabulate.cpp input.cpp elonO2.cpp
	$(CC) $(CFLAGS) -o $@ test_scatter_speed.cpp gnuplot_i.c histogram.cpp $(OPTIM) $(INC) $(LIB) 
test_flux: test_flux.cpp random.c random.h
	$(CC) $(CFLAGS) $(OPTIM) -o ../$@ test_flux.cpp random.c histogram.cpp
rot_test: rot_test.cpp random.c random.h
	$(CC) $(CFLAGS) $(OPTIM) -o ../$@ rot_test.cpp random.c
maxwell_test: maxwell_test.cpp random.c random.h
	$(CC) $(CFLAGS) $(OPTIM) -o ../$@ maxwell_test.cpp random.c gnuplot_i.c histogram.cpp
gnuplot_i.o: gnuplot_i.h
histogram.o: histogram.hpp
