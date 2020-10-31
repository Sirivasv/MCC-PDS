all:
	g++ -g -O2 -I/usr/include/ ReadMicWavs.cpp -o ReadMicWavs -ljack -lsndfile
	g++ MDOAMUSIC.cpp -o out_MDOAMUSIC -ljack -lfftw3 -lm -I/usr/include/eigen3
testEigen:
	g++ -I/usr/include/eigen3 eigen_examples.cpp -o eigen_examples

