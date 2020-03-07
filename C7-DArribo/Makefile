all:
	g++ -g -O2 -I/usr/include/ ReadMicWavs.cpp -o ReadMicWavs -ljack -lsndfile
	gcc darribo.c -o darribo -ljack -lfftw3 -lm
testEigen:
	g++ -I/usr/include/eigen3 eigen_examples.cpp -o eigen_examples

