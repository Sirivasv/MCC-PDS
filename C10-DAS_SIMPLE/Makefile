all:
	g++ -g -O2 -I/usr/include/ ReadMicWavs.cpp -o ReadMicWavs -ljack -lsndfile
	g++ DAS_SIMPLE.cpp -o out_DAS_SIMPLE -lsndfile -ljack -lfftw3 -lm -I/usr/include/eigen3
testEigen:
	g++ -I/usr/include/eigen3 eigen_examples.cpp -o eigen_examples

