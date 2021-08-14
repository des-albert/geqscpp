CC = g++
LD = g++

CFLAGS = -g 

SOURCES = geq.cpp bndmat.cpp eqsil.cpp flux.cpp intpol.cpp splnco.cpp startt.cpp compar.cpp xcur.cpp gauss.cpp saddle.cpp \
	curf.cpp topol.cpp plotit.cpp

OBJECTS = $(SOURCES:.cpp=.o)

.SUFFIXES: .cpp .o 
%.o: %.cpp
	$(CC) $(CFLAGS) -c $<


gep: $(OBJECTS)
	$(LD) -o geq $(OBJECTS) -L/usr/local/dislin -ldiscpp
all: gep
clean: 
	rm -f *.o geq
