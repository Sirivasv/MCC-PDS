/* Some examples of how to use Eigen library */

/* Compile with: g++ -I/usr/include/eigen3 eigen_examples.cpp -o eigen_examples */

#include <stdio.h>		//for printf
#include <iostream>		//for std::cout, as a replacemente for printf

// Eigen include
#include <Eigen/Eigen>


int main() {
	int i,j;
	
	/*
	 * The following defines A as a 2x2 matrix, filled with doubles.
	 * That's why it ends with 'd'. If you want to use floats,
	 * use MatrixXf.
	 */
	
	Eigen::MatrixXd A(2,2);
	A(0,0) = 1.0; A(0,1) = 2.0; A(1,0) = 2.0; A(1,1) = 3.0;
	
	printf("--- Matrix A:\n");
	for(i=0; i < 2; i++){
		for(j=0; j < 2; j++){
			printf("\t %0.2f",A(i,j));
		}
		printf("\n");
	}
	printf("\n");
	/*
	 * If iostream is used instead of stdio.h, the following line can
	 * replace the last two for's that print the A matrix:
	 */
	//std::cout << "The matrix A is:\n" << A << std::endl << std::endl;

	
	/*
	 * EigenValues and the components of the EigenVectors are usually 
	 * returned as std::complex type.
	 */
	
	Eigen::EigenSolver<Eigen::MatrixXd> es(A);
	printf("--- Eigenvalues of A:\n");
	for(i=0; i < 2; i++){
		printf("\t%d: \t%f\n", i+1,es.eigenvalues()[i].real());
	}
	printf("\n");
	/*
	 * If iostream is used instead of stdio.h, the following line can
	 * replace the last two for's that print the EigenValues list:
	 */
	//std::cout << "The eigenvalues of A are:\n" << es.eigenvalues().real() << std::endl << std::endl;
	
	printf("--- Eigenvectors of A:\n");
	for(i=0; i < 2; i++){
		printf("\t%d:", i+1);
		for(j=0; j < 2; j++){
			printf("\t %f",es.eigenvectors().col(i)[j].real());
		}
		printf("\n");
	}
	printf("\n");
	/*
	 * If iostream is used instead of stdio.h, the following line can
	 * replace the last two for's that print the EigenVectors:
	 */
	//std::cout << "The eigenvectors of A are (one vector per column):\n" << es.eigenvectors().real() << std::endl << std::endl;
	
	
	/*
	 * For simplicity, from this point on, std::cout will be used to print.
	 */
	
	Eigen::MatrixXd B(2,2);

	B.setZero();
	std::cout << "--- B.setZero():\n" << B << std::endl << std::endl;

	B.setOnes();
	std::cout << "--- B.setOnes():\n" << B << std::endl << std::endl;
	
	B.setRandom();
	std::cout << "--- B.setRandom():\n" << B << std::endl << std::endl;
	
	B.setConstant(2.5);
	std::cout << "--- B.setConstant(2.5):\n" << B << std::endl << std::endl;
	
	std::cout << "--- A+B:\n" << A+B << std::endl << std::endl;
	std::cout << "--- A-B:\n" << A-B << std::endl << std::endl;
	std::cout << "--- A.cwiseProduct(B) == Matlab -> A.*B :\n" << A.cwiseProduct(B) << std::endl << std::endl;
	std::cout << "--- A.cwiseQuotient(B) == Matlab -> A./B :\n" << A.cwiseQuotient(B) << std::endl << std::endl;
	std::cout << "--- A*B:\n" << A*B << std::endl << std::endl;
	std::cout << "--- A.row(0)*B == Matlab -> A(1,:)*B :\n" << A.row(0)*B << std::endl << std::endl;
	std::cout << "--- A*B.col(1) == Matlab -> A*B(:,2) :\n" << A*B.col(1) << std::endl << std::endl;
	std::cout << "--- A.transpose() == Matlab -> A' :\n" << A.transpose() << std::endl << std::endl;
	std::cout << "--- A.inverse() == Matlab -> inv(A) :\n" << A.inverse() << std::endl << std::endl;
	
	/*
	 * Eigen uses the standard implementation of complex number:
	 * std::complex.
	 * 
	 * In this examples, the last 'd' in MatrixXcd means that C is a 
	 * matrix that is full of complex numbers that use doubles for its
	 * real and imaginary parts and thus it needs: std::complex<double>.
	 * 
	 * To use floats instead: use MatrixXcf and std::complex<float>.
	 */ 
	Eigen::MatrixXcd C(2,2);
	C(0,0) = std::complex<double>(1.0,0.5);
	C(0,1) = std::complex<double>(2.0,1.6);
	C(1,0) = std::complex<double>(3.0,2.7);
	C(1,1) = std::complex<double>(4.0,3.8);
	std::cout << "--- Complex Number Matrix C:\n" << C << std::endl << std::endl;
	
	C(0,0).real(5.0);
	std::cout << "--- C(0,0).real(5.0):\n" << C << std::endl << std::endl;
	
	Eigen::ComplexEigenSolver<Eigen::MatrixXcd> es_c(C);
	std::cout << "--- The eigenvalues of C are:\n" << es_c.eigenvalues() << std::endl << std::endl;
	std::cout << "--- The eigenvectors of C are (one vector per column):\n" << es_c.eigenvectors() << std::endl << std::endl;
	
	return 0;
}
