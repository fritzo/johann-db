
current_target: test

CXX = g++
CC = $(CXX)
CPPFLAGS = -I ../lib
CXXFLAGS = -fno-exceptions -fno-rtti \
	   -Wall -Wextra -Winit-self -pedantic -Wno-unused-parameter -ggdb \
	   -O3 -fopenmp -pipe -march=native -ffast-math \
	   -funsafe-math-optimizations -ffinite-math-only -fno-trapping-math \
	   -mrecip -mfpmath=sse -fvect-cost-model
CFLAGS = $(CXXFLAGS)
LDFLAGS = -lm -lpthread

jdb.o: ../lib/jdb.cpp ../lib/jdb.h
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c $< -o $@

jdb_to_csv: jdb_to_csv.o jdb.o

test: jdb_to_csv
	./jdb_to_csv ~/johann/data/default.jdb ./test

all: jdb_to_csv

clean:
	rm -f *.o jdb_to_csv *.csv

