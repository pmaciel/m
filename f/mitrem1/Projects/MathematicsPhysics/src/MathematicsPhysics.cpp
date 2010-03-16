//---------------------------------------------------------------------------

#include "MathematicsPhysics.h"
#include <algorithm>
#include <iostream>
#include <fstream>
#include <string>

//---------------------------------------------------------------------------

void errorSingularMatrix();
void errorZero(const std::string &cppfile, const std::string &function);

//---------------------------------------------------------------------------
// Faculty
double fac(double N)
{
	double facN = 1;
	for (double i=2; i<=N; i++)
		facN *= i;
	return facN;
}
//---------------------------------------------------------------------------
// Kronecker symbol
unsigned kroneck(unsigned i, unsigned j)
{
	unsigned k = 0;
	if (i==j) k = 1;
	return k;
}
//---------------------------------------------------------------------------
// Exponential integral
double expInt(double x)
{
	double integral = 0.;
	double factor;
	double term;
	if (x <= 0) errorZero("Mathematics.cpp","expInt");
	else if (x > 7) {
		double k = 1;
		integral = 1.;
		do {
			factor = 1.;
			for (unsigned p=0; p<k; p++)
				factor *= -1./x;
			term = fac(k)*factor;
			integral += term;
			k++;
		} while ((fabs(term/integral) > 1e-2) && (k<7)); // Convergence will be reached before k=6 if x > 7
		integral *= exp(-x)/x;
	}
	else {
		integral = -eulergamma_CONST-log(x);
		double k = 1;
		do {
			factor = -x;
			for (unsigned p=1; p<k; p++)
				factor *= (-x);
			term = factor/(k*fac(k));
			integral -= term;
			k++;
		} while ((fabs(term/integral) > 1e-2) && (k<27)); // Convergence will be reached before k=26 if x < 7
	}
	return integral;
}
//---------------------------------------------------------------------------
// Gauss solver with pivoting
void solveGauss(double** A, double* B, unsigned N)
{
	double C;
	for (unsigned m=0; m<N-1; m++) {
		// Put row with highest diagonal element on top
		C = A[m][m];
		for (unsigned n=m+1; n<N; n++) {
			if (fabs(A[n][m]) > fabs(C)) {
				for (unsigned p=m; p<N; p++) {
					C = A[m][p];
					A[m][p] = A[n][p];
					A[n][p] = C;
				}
				C = B[m];
				B[m] = B[n];
				B[n] = C;
				C = A[m][m];
			}
		}

		// Check if diagonal element is (close to) zero
		if (fabs(C) < 1.0e-32) errorSingularMatrix();

		// Normalize row m
		for (unsigned n=m+1; n<N; n++)
			A[m][n] /= C;
		B[m] /= C;

		// Subtract row m from subsequent rows 
		for (unsigned n=m+1; n<N; n++) {
			C = A[n][m];
			for (unsigned p=m+1; p<N; p++)
				A[n][p] -= C*A[m][p];
			B[n] -= C*B[m];
		}
		
		//cout << '-';
	}

	// Solve by back substitution
	B[N-1] /= A[N-1][N-1];
	for (unsigned p=0; p<N-1; p++) {
		unsigned m = N-p-2;
		for (unsigned n=m+1; n<N; n++)
			B[m] -= A[m][n]*B[n];
	}
}
//---------------------------------------------------------------------------


//--- ERROR MESSAGES --------------------------------------------------------
void errorZero(const std::string &cppfile, const std::string &function)
{
	std::cout << "ERROR IN " << cppfile << ".\nA PARAMETER IN " << function << " IS <= 0." << std::endl;
	system("pause");
	exit(1);
}
//---------------------------------------------------------------------------
void errorSingularMatrix()
{
	std::cout << "Mathematics.cpp\".\
			\nTHE MATRIX IS SINGULAR." << std::endl;
	system("pause");
	exit(1);
}
//---------------------------------------------------------------------------

