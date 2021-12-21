all: par trayectorias poblaciones grav

par:
	g++ -std=c++11 -O3  SOURCE/par.cpp SOURCE/random.cpp -o BIN/par.bin
trayectorias:
	g++ -std=c++11 -O3  SOURCE/trayectorias.cpp SOURCE/random.cpp -o BIN/trayectorias.bin
poblaciones:
	g++ -std=c++11 -O3  SOURCE/poblaciones.cpp SOURCE/random.cpp -o BIN/poblaciones.bin
grav:
	g++ -std=c++11 -O3  SOURCE/interaccion_grav.cpp SOURCE/random.cpp -o BIN/grav.bin

clean:
	rm -rf BIN/*.bin
