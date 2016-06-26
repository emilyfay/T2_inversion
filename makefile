#makefile for CPMG inversion

CC = g++
CFLAGS = -g -Wall -c
LFLAGS = -Wall -g
OBJS = CPMG.o T2_inversion.o NNLS.o base_funcs.o
INCLUDES = -I/usr/include/

default: CPMG

#to create executable CPMG
CPMG: $(OBJS)
	$(CC) $(LFLAGS) $(INCLUDES) -o CPMG $(OBJS)  -O2 -DARMA_DONT_USE_WRAPPER -framework Accelerate
	
#to create object file for CPMG
CPMG.o: CPMG.cpp T2_inversion.h
	$(CC) $(CFLAGS) $(INCLUDES)  CPMG.cpp -O2 -DARMA_DONT_USE_WRAPPER
#to create object file for T2_inversion
T2_inversion.o: T2_inversion.h T2_inversion.cpp NNLS.h base_funcs.h
	$(CC) $(CFLAGS) $(INCLUDES)  T2_inversion.cpp -O2 -DARMA_DONT_USE_WRAPPER
#to create object file for NNLS
NNLS.o: NNLS.h NNLS.cpp base_funcs.h
	$(CC) $(CFLAGS) $(INCLUDES) NNLS.cpp  -O2 -DARMA_DONT_USE_WRAPPER #-framework Accelerate
#to create object file for base_funcs
base_funcs.o: base_funcs.h base_funcs.cpp
	$(CC) $(CFLAGS) base_funcs.cpp
	
clean:
	\rm *.o *~ CPMG
	