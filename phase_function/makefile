#
mie: main.o mie.o
	g++ -o mie main.o mie.o
main.o: main.cpp mie.h
	g++ -c main.cpp
mie.o: mie.cpp mie.h
	g++ -c mie.cpp
clean:
	rm -f mie main.o mie.o
