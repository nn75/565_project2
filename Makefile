CFLAGS = -O3 -Wall -pedantic -std=c++14
CC = g++
LIB = -lpthread

rainfall: rainfall_seq.cpp rainfall_pt.cpp RainfallSimulator.o
	$(CC) $(CFLAGS) -o rainfall_seq rainfall_seq.cpp RainfallSimulator.o
	$(CC) $(CFLAGS) -o rainfall_pt rainfall_pt.cpp $(LIB)

RainfallSimulator.o: RainfallSimulator.cpp RainfallSimulator.h
	$(CC) $(CFLAGS) -c RainfallSimulator.cpp
clean:
	rm -f rainfall_seq rainfall_pt *.o *.c~ *.h~ *.cpp~
