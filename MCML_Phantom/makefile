#
mc: main.o io.o transport.o fiber.o
	nvcc -o mc main.o io.o transport.o fiber.o
main.o: main.cu header.h
	nvcc -c main.cu
io.o: io.cu header.h
	nvcc -c io.cu
transport.o: transport.cu header.h
	nvcc -c transport.cu
fiber.o: fiber.cu header.h
	nvcc -c fiber.cu
clean:
	rm -f mc main.o io.o transport.o fiber.o
