//---------------------------------------------------------------------------
#include "SolverLinear.h"
#include <math.h>
#include <algorithm>
#include <iostream>
#include <fstream>
#include "TypeDefs.h"

#include "mkl.h"

//---------------------------------------------------------------------------


//---------------------------------------------------------------------------
SolverLinear::SolverLinear (const Element* const* elements, unsigned nElements, unsigned nNodes, unsigned nVariables)
{
	nRows = nNodes*nVariables;
	
	std::vector< std::vector<unsigned> >::iterator r;
	std::vector<unsigned>::iterator c;

	// Determine matrix structure (store the nodes to which each node is connected)
	matStruct.resize(nNodes);
	bool alreadyIncluded;
	for (unsigned e=0; e<nElements; e++) 
	{
		const Element* element = elements[e];
		unsigned nElementNodes = element->getNNodes();
		IndexList nodes = element->getNodes();
		for (unsigned m=0; m<nElementNodes; m++) 
		{
			for (unsigned n=0; n<nElementNodes; n++) 
			{
				alreadyIncluded = false;
				for (c = matStruct[nodes[m]].begin(); c != matStruct[nodes[m]].end(); c++) 
				{
					if (*c == nodes[n]) 
					{
						alreadyIncluded = true;
						break;
					}
				}
				if (!alreadyIncluded)
				matStruct[nodes[m]].push_back(nodes[n]);
			}
		}
	}

	// Sort the column indices
	for (r = matStruct.begin(); r != matStruct.end(); r++) 
	{
		sort(r->begin(),r->end());
	}
}
//---------------------------------------------------------------------------


//--- BAND LU ---------------------------------------------------------------
SolverLinear_BandLU::SolverLinear_BandLU (const Element* const* elements, unsigned nElements, unsigned nNodes, unsigned nVariables) 
	: SolverLinear (elements, nElements, nNodes, nVariables)
{
	std::vector< std::vector<unsigned> >::iterator r;
	std::vector<unsigned>::iterator c;

	// Determine number of lower and upper co-diagonals
	ld = 0;
	ud = 0;
	unsigned row = 0;
	unsigned min, max;
	for (r = matStruct.begin(); r != matStruct.end(); r++) 
	{
		min = *r->begin();
		max = *r->begin();
		for (c = r->begin(); c != r->end(); c++) 
		{
			if (*c < min) min = *c;
			if (*c > max) max = *c;
		}
		if (max-row > ld) ld = max - row;
		if (row-min > ud) ud = row - min;
		row++;
	}
	LD = (ld+1)*nVariables-1;
	UD = (ud+1)*nVariables-1;
	nColumns = LD+1+UD;

	//ofstream Band("Bandwidth.out");
	//Band << "LD = " << LD << "\tUD = " << UD << endl;
	//Band << "Bandwidth = " << double(LD+1+UD)/((nIons+1)*nNodes)*100 << endl;

	// Build matrix A
	A = new double*[nRows];
	for (unsigned row=0; row<nRows; row++)
	{
		A[row] = new double[nColumns];
	}

}
//---------------------------------------------------------------------------
SolverLinear_BandLU::~SolverLinear_BandLU () 
{
	for (unsigned row=0; row<nRows; row++)
	{
		delete[] A[row];
	}
	delete[] A;
}
//---------------------------------------------------------------------------
void SolverLinear_BandLU::init()
{
	for (unsigned row=0; row<nRows; row++)
	{
		for (unsigned col=0; col<nColumns; col++)
		{
			A[row][col] = 0.;
		}
	}
}
//---------------------------------------------------------------------------
void SolverLinear_BandLU::add (unsigned row, unsigned col, double value)
{
	A[row][LD+col-row] += value;
}
//---------------------------------------------------------------------------
void SolverLinear_BandLU::multiplyWithVector (double* X, double* Result)
{
	for (unsigned row=0; row<nRows; row++)
	{
		double sum = 0.;
		for (unsigned col=row-std::min(row,LD); col<std::min(nRows,UD+1+row); col++)
		{
			sum += A[row][LD+col-row]*X[col];
		}
		//cout << sum << endl;
		Result[row] = sum;
	}
	
	//string f = "Test.out";
	//Print(f);
	// You may replace var(m,i) by eq(m,i) to be more logical, but it will run slower.
	//for (unsigned m=0; m<nNodes; m++) 
	//{
	//	for (unsigned i=0; i<nVariables; i++) 
	//	{
	//		unsigned varmi = var(m,i);
	//		double sum = 0.;
	//		for (unsigned n=m-std::min(m,ld); n<std::min(nNodes,ud+1+m); n++) 
	//		{
	//			for (unsigned j=0; j<nVariables; j++) 
	//			{
	//				//cout << A[var(m,i)][LD+var(n,j)-var(m,i)] << '\t' << X[var(n,j)] << endl;
	//				sum += A[varmi][LD+var(n,j)-varmi]*X[var(n,j)];
	//			}
	//		}
	//		//cout << sum << endl;
	//		Result[varmi] = sum;
	//	}
	//}
}
//---------------------------------------------------------------------------
void SolverLinear_BandLU::residu (double* X, double* B)
{
	for (unsigned row=0; row<nRows; row++)
	{
		double sum = 0.;
		for (unsigned col=row-std::min(row,LD); col<std::min(nRows,UD+1+row); col++)
		{
			sum += A[row][LD+col-row]*X[col];
		}
		//cout << sum << endl;
		B[row] -= sum;
	}
	//string f = "Test.out";
	//Print(f);
	// You may replace var(m,i) by eq(m,i) to be more logical, but it will run slower.
	//for (unsigned m=0; m<nNodes; m++) 
	//{
	//	for (unsigned i=0; i<nVariables; i++) 
	//	{
	//		unsigned varmi = var(m,i);
	//		double sum = 0.;
	//		for (unsigned n=m-std::min(m,ld); n<std::min(nNodes,ud+1+m); n++) 
	//		{
	//			for (unsigned j=0; j<nVariables; j++) 
	//			{
	//				//cout << A[var(m,i)][LD+var(n,j)-var(m,i)] << '\t' << X[var(n,j)] << endl;
	//				sum += A[varmi][LD+var(n,j)-varmi]*X[var(n,j)];
	//			}
	//		}
	//		//cout << sum << endl;
	//		B[varmi] -= sum;
	//	}
	//}
}
//---------------------------------------------------------------------------
double SolverLinear_BandLU::getResidual (double* X, double* B, unsigned int row)
{
	double residual = 0.;

	for (unsigned col=row-std::min(row,LD); col<std::min(nRows,UD+1+row); col++)
	{
		residual += A[row][LD+col-row]*X[col];
	}
	residual -= B[row];
	return residual;
}
//---------------------------------------------------------------------------
void SolverLinear_BandLU::imposeVar (unsigned row)
{
	// Here, you may NOT replace var(m,i) by eq(m,i)!!!
	for (unsigned col=row-std::min(row,LD); col<std::min(nRows,UD+1+row); col++)
	{
		A[row][LD+col-row] = 0.;
	}
	A[row][LD] = 1.;
	
	//unsigned varmi = var(m,i);
	//// Here, you may NOT replace var(m,i) by eq(m,i)!!!
	//for (unsigned n=m-std::min(m,ld); n<std::min(nNodes,ud+1+m); n++) 
	//{
	//	for (unsigned j=0; j<nVariables; j++) 
	//	{
	//		A[varmi][LD+var(n,j)-varmi] = 0.;
	//	}
	//}
	//A[varmi][LD] = 1.;
}
//---------------------------------------------------------------------------
void SolverLinear_BandLU::clearRow (unsigned row)
{
	for (unsigned col=row-std::min(row,LD); col<std::min(nRows,UD+1+row); col++)
	{
		A[row][LD+col-row] = 0.;
	}
}
//---------------------------------------------------------------------------
void SolverLinear_BandLU::linearCombination (unsigned toAddTo, unsigned toAdd, double factor)
{
	for (unsigned col=toAddTo-std::min(toAddTo,LD); 
				col<std::min(nRows,UD+1+toAddTo); col++)
	{
		A[toAddTo][LD+col-toAddTo] += factor*A[toAdd][LD+col-toAdd];
	}
}
//---------------------------------------------------------------------------
void SolverLinear_BandLU::LUDecomp (double* B)

/*=========================================================================*
 *                                                                         *
 *  The function LUDecomp factors a condensed banded matrix A.             *
 *  LUDecomp uses the Gauss algorithm without column pivot search.         *
 *                                                                         *
 *=========================================================================*/
{
	unsigned kend, kjend, jk, jm;

	// Loop over all rows.
	for (unsigned i=0; i<nRows-1; i++) {
		kend  = std::min (LD+1, nRows-i);
		kjend = std::min (UD+1, nRows-i);

		// Check if matrix is singular.
		if (fabs(A[i][LD]) < 1e-32) errorSingularMatrix();

		// Loop over all rows below row i.
		for (unsigned k=1; k<kend; k++) {
			A[k+i][LD-k] /= A[i][LD];
			B[k+i] -= A[k+i][LD-k]*B[i];

			// Loop over columns.
			for (unsigned j=1; j!=kjend; j++) {
				jk = j+LD-k;
				jm = j+LD;
				A[k+i][jk] -= A[k+i][LD-k]*A[i][jm];
			}

		}

	}
}
//---------------------------------------------------------------------------
void SolverLinear_BandLU::solve (double* B)

/*=========================================================================*
 *                                                                         *
 *  The function Solve solves a linear banded system: A * X = B.           *
 *  - A is a nonsingular size x size matrix in                             *
 *    condensed form, i.e. represented in a                                *
 *    size x (ld+ud+1)*nVariables matrix.                                  *
 *  - B denotes the right hand side of the system.                         *
 *  - X is the solution.                                                   *
 *                                                                         *
 *=========================================================================*
 *                                                                         *
 *   Applications:                                                         *
 *   =============                                                         *
 *      Solve linear systems with nonsingular banded system matrices.      *
 *      Particularly useful for large sparse and banded and diagonally     *
 *      dominant matrices.                                                 *
 *                                                                         *
 *=========================================================================*
 *                                                                         *
 *   Input parameters:                                                     *
 *   ================                                                      *
 *      (nIons+1)*nNodes        Dimension of A,                            *
 *                              size of B ((nIons+1)*nNodes > 2)           *
 *      ld      Number of lower co-diagonals (ld >= 0)                     *
 *      ud      Number of upper co-diagonals (ud >= 0)                     *
 *      A       Matrix of the system in comdensed form.                    *
 *              Each row has length at least ld + 1 + ud + std::min(ld,ud)      *  (+ std::min(ld,ud)? why?)
 *              where the columns 0, .., ld-1 denote the lower             *
 *              co-diagonals, column ld stores the diagonal and the        *
 *              columns ld+1, .., ld+ud contain the upper                  *
 *              co-diagonals.                                              *
 *              If AMat is the original uncondensed band matrix:           *
 *                 AMat[i][k] = A[i][ld+k-i],                              *
 *                 for k,i inside the band                                 *
 *      B       Right hand side                                            *
 *                                                                         *
 *   Output parameters:                                                    *
 *   ==================                                                    *
 *      A       LU factorization in condensed form                         *
 *      B       solution vector for the system                             *
 *                                                                         *
 *=========================================================================*/
{
	// Check if matrix is not already in upper triangular form.
	if (ld > 0) {
		LUDecomp(B);
	}

	unsigned kend;

	// Back substitution
	for (int i=nRows-1; i>=0; i--) {
		kend = std::min (UD+1, nRows-i);
		for (unsigned k=1; k<kend; k++) {
			B[i] -= A[i][k+LD] * B[i+k];
		}
		B[i] /= A[i][LD];
	}

}
//---------------------------------------------------------------------------
void SolverLinear_BandLU::print(const std::string &filename)
{
	std::ofstream CHECK;
	CHECK.open(filename.c_str(), std::ofstream::out | std::ofstream::app);
	CHECK.precision(12);
	/*CHECK << "A form : " << std::endl;
	for (unsigned row=0; row<nRows; row++) {
		for (unsigned col=0; col<nColumns; col++) {
			CHECK << A[row][col] << '\t';
		}
		CHECK << std::endl;
	}*/
	CHECK << "A : " << std::endl;
	for (unsigned row=0; row<nRows; row++) {
		for (unsigned col=0; col<nRows; col++) {
			if ((LD+col-row < nColumns) && (LD+col-row >= 0))
			{
				CHECK << A[row][LD+col-row] << '\t';
			}
			else
			{
				CHECK << "0" << "\t";
			}
		}
		CHECK << std::endl;
	}
	CHECK.close();

	//std::ofstream CHECK;
	//CHECK.open(filename.c_str(), std::ofstream::out | std::ofstream::app);
	//CHECK.precision(12);
	//CHECK << "A form : " << std::endl;
	//for (unsigned m=0; m<nNodes; m++) {
	//	for (unsigned i=0; i<nVariables; i++) {
	//		for (unsigned n=0; n<nNodes; n++) {
	//			for (unsigned j=0; j<nVariables; j++) {
	//				if ((n<=m+ud) && (m<=n+ld)) {
	//					CHECK << A[var(m,i)][LD+var(n,j)-var(m,i)] << "\t";
	//					//cout << A[var(m,i)][LD+var(n,j)-var(m,i)] << "\t";
	//				}
	//				else {
	//					CHECK << "0" << "\t";
	//					//cout << "0" << "\t";
	//				}
	//			}
	//		}
	//		CHECK << std::endl;
	//		//cout << endl;
	//	}
	//}
	//CHECK << std::endl;
	////cout << endl;
	//CHECK.close();

	/*CHECK << "P form : " << endl;
	for (unsigned m=0; m<(SNLD->mitrem->data->nIons+1)*SNLD->reactor->nNodes; m++) {
		for (unsigned n=0; n<LD+1+UD; n++) {
			//if ( (n+m>=LD) && (n+m<=LD+2+UD) )
				CHECK << A[m][n] << "\t";
			//else
				//output << "0" << "\t";
		}
		CHECK << endl;
	}
	CHECK << endl;*/

	/*cout << "A form : " << endl;
	for (unsigned m=0; m<SNLD->reactor->nNodes; m++) {
		for (unsigned i=0; i<SNLD->mitrem->data->nIons+1; i++) {
			for (unsigned n=0; n<SNLD->reactor->nNodes; n++) {
				for (unsigned j=0; j<SNLD->mitrem->data->nIons+1; j++) {
				if ((n<=m+ud) && (m<=n+ld))
					cout << A[SNLD->var(m,i)][LD+SNLD->var(n,j)-SNLD->var(m,i)] << "\t";
				else
					cout << "0" << "\t";
				}
			}
			cout << endl;
		}
	}
	cout << endl;
	stop();*/
}
//---------------------------------------------------------------------------
void SolverLinear_BandLU::errorSingularMatrix()
{
	std::cout << "BandSolver.cpp\".\
			\nTHE MATRIX IS SINGULAR." << std::endl;
	system("pause");
	exit(1);
}
//---------------------------------------------------------------------------


//--- GMRES -----------------------------------------------------------------
SolverLinear_GMRES::SolverLinear_GMRES (const Element* const* elements, unsigned nElements, unsigned nNodes, unsigned nVariables) 
	: SolverLinear (elements, nElements, nNodes, nVariables)
{
	this->nVariables = nVariables;
	nRows = nNodes*nVariables;
	
	std::vector< std::vector<unsigned> >::iterator r;
	std::vector<unsigned>::iterator c;

	// Determine the number of non-zero's
	N = 0;
	for (r = matStruct.begin(); r != matStruct.end(); r++) 
	{
		N += r->size();
	}

	// Build matrix A
	Rows = new int[nRows+1];
	Columns = new int[N*nVariables*nVariables];
	A = new double[N*nVariables*nVariables];

	unsigned rowNode, colNode, rowIndex, colIndex;

	Rows[0] = 1;
	   
	rowNode = 0;
	for (r = matStruct.begin(); r != matStruct.end(); r++) 
	{
		for (unsigned rowVar=0; rowVar<nVariables; rowVar++) 
		{
			rowIndex = rowNode*nVariables+rowVar;
			Rows[rowIndex+1] = Rows[rowIndex] + r->size()*nVariables;
			colNode = 0;
			for (c = r->begin(); c != r->end(); c++) 
			{
				for (unsigned colVar=0; colVar<nVariables; colVar++) 
				{
					colIndex = Rows[rowIndex]+colNode*nVariables+colVar-1;
					Columns[colIndex] = (*c)*nVariables+colVar+1;
				}
				colNode++;
			}
		}
		rowNode++;
	}

	c__1 = 1;

	for (unsigned i=0;i<64;++i)
	{
		pt_[i] = 0;
	}

}
//---------------------------------------------------------------------------
SolverLinear_GMRES::~SolverLinear_GMRES () 
{
	delete[] Rows;
	delete[] Columns;
	delete[] A;

}
//---------------------------------------------------------------------------
void SolverLinear_GMRES::init()
{
	for (unsigned i=0; i<N*nVariables*nVariables; i++)
	{
		A[i] = 0;
	}
}
//---------------------------------------------------------------------------
unsigned SolverLinear_GMRES::getIndex(unsigned row, unsigned col)
{
	unsigned index = Rows[row]-1;
	while ((Columns[index] != col+1) && (index < Rows[row+1]-1))
		index++;
	return index;
}
//---------------------------------------------------------------------------
void SolverLinear_GMRES::add (unsigned row, unsigned col, double value)
{
	A[getIndex(row,col)] += value;
}
//---------------------------------------------------------------------------
void SolverLinear_GMRES::multiplyWithVector (double* X, double* Result)
{
	for (unsigned r=0; r<nRows; r++) 
	{
		double sum = 0;
		for (unsigned c=Rows[r]-1; c<Rows[r+1]-1; c++) 
		{
			sum += A[c]*X[Columns[c]-1];
		}
		Result[r] = sum;
	}
}
//---------------------------------------------------------------------------
void SolverLinear_GMRES::residu (double* X, double* B)
{
	for (unsigned r=0; r<nRows; r++) 
	{
		double sum = 0;
		for (unsigned c=Rows[r]-1; c<Rows[r+1]-1; c++) 
		{
			sum += A[c]*X[Columns[c]-1];
		}
		B[r] -= sum;
	}
}
//---------------------------------------------------------------------------
double SolverLinear_GMRES::getResidual (double* X, double* B, unsigned int row)
{
	double residual = 0.;
	for (unsigned c=Rows[row]-1; c<Rows[row+1]-1; c++) 
	{
		residual += A[c]*X[Columns[c]-1];
	}
	residual -= B[row];
	return residual;
}
//---------------------------------------------------------------------------
void SolverLinear_GMRES::imposeVar (unsigned row)
{
	for (unsigned c=Rows[row]-1; c<Rows[row+1]-1; c++) 
	{
		A[c] = 0.;
	}
	A[getIndex(row,row)] = 1.;
}
//---------------------------------------------------------------------------
void SolverLinear_GMRES::clearRow (unsigned row)
{
	for (unsigned c=Rows[row]-1; c<Rows[row+1]-1; c++) 
	{
		A[c] = 0.;
	}
}
//---------------------------------------------------------------------------
void SolverLinear_GMRES::linearCombination (unsigned toAddTo, unsigned toAdd, double factor)
{
	unsigned c2 = Rows[toAdd]-1;
	for (unsigned c=Rows[toAddTo]-1; c<Rows[toAddTo+1]-1; c++) 
	{
		A[c] += factor*A[c2];
		c2++;
	}
}
//---------------------------------------------------------------------------
void SolverLinear_GMRES::solve(double* B)
{
	solveWithPardiso(B);
}
//---------------------------------------------------------------------------
void SolverLinear_GMRES::solveWithGMRES (double* B)
{
	int ierr;
	int n = nRows;
	int iwk = N*nVariables*nVariables;
	int im = 50;
	double eps = 1e-5;
	int maxits = 50;
	int iout = 0;
	int lfil = 3;

	double *alu = new double[iwk];
	int *jlu = new int[iwk];
	int *levs = new int[iwk];
	for (register int i=0; i<iwk; i++) {
		alu[i] = 0.;
		jlu[i] = 0;
	}

	int *ju = new int[n+1];
	for (register int i=0; i<n; i++) {
		ju[i] = 0;
	}

	double *w = new double[n+1];
	int *jw = new int[n*3];
	for (register int i=0; i<n; i++) {
		w[i] = 0.;
		jw[i] = 0;
		jw[n+i] = 0;
		jw[2*n+i] = 0;
	}

	int newiwk = iluk(&n,A,(int*)Columns,(int*)Rows,&lfil,
	alu,jlu,ju,levs,&iwk,w,jw,&ierr);
	   
	delete[] w;
	delete[] jw;
	delete[] levs;

	double *vv = new double[n*(im+1)];
	double *dblRHS = new double[n];
	double *dblRES = new double[n];

	for (unsigned index=0; index<n; index++) {
		dblRHS[index] = B[index];
		dblRES[index] = 0;
	}

	pgmres(&n,&im,dblRHS,dblRES,vv,&eps,&maxits,&iout,
				A,(int*)Columns,(int*)Rows,alu,jlu,ju,&ierr);

	for (unsigned index=0; index<n; index++) {
		B[index] = dblRES[index];
	}

	delete [] alu;
	delete [] jlu;
	delete [] ju;
	delete [] vv;
	delete [] dblRHS;
	delete [] dblRES;
}
//---------------------------------------------------------------------------
void SolverLinear_GMRES::solveWithPardiso (double* B)
{
	int n = nRows;

	int maxfct = 1;
	int mnum = 1;
	//int mtype = bridge_->linsolvercontrol_.linSolver_.pardiso_mtype_;
	//int phase = bridge_->linsolvercontrol_.linSolver_.pardiso_phase_;
	int mtype = 11;

	//static int first = 0;
	int phase = 13;
	//if (first<10)
	//{
	//	phase = 13;
	//	first++;
	//}

	int* perm = 0;

	int iparm[64];
	for (unsigned i=0;i<64;++i)
	{
		iparm[i] = 0;
	}
	//iparm[0] = 1;
	iparm[5] = 1;

	int error = 0;

	int one = 1;
	int msglevel = 0;

	//double *dblRHS = new double[SNLD->reactor->nNodes*(SNLD->mitrem->data->nIons+1)];
	double *dblRES = new double[n];

	for (unsigned index=0; index<n; index++) {
		//dblRHS[index] = B[index];
		dblRES[index] = 0;
	}

	/*std::ofstream matrix("Pardiso.out");
	matrix << "IA\n";
	for (unsigned r=0; r<nRows+1; r++)
	{
		matrix << r << '\t' << Rows[r] << std::endl;
	}
	matrix << "JA\n";
	for (unsigned c=0; c<N*nVariables*nVariables; c++)
	{
		matrix << c << '\t' << Columns[c] << std::endl;
	}
	matrix << "A\n";
	for (unsigned a=0; a<N*nVariables*nVariables; a++)
	{
		matrix << a << '\t' << A[a] << std::endl;
	}
	matrix << "B\n";
	for (unsigned b=0; b<n; b++)
	{
		matrix << b << '\t' << B[b] << std::endl;
	}
	matrix << std::endl;*/

	std::cout << "Running Pardiso..." << std::endl;
	PARDISO(pt_,&maxfct,&mnum,&mtype,&phase,&n,A,Rows,Columns,perm,&one,iparm,&msglevel,B,dblRES,&error);
	std::cout << "Pardiso has finished..." << std::endl;

	/*matrix << "Result\n";
	for (unsigned b=0; b<n; b++)
	{
		matrix << b << '\t' << B[b] << std::endl;
	}
	matrix << std::endl;
	exit(1);*/

	if (error != 0)
		std::cout << "Error for Pardiso = " << error << std::endl;

	for (unsigned index=0; index<n; index++) {
		B[index] = dblRES[index];
	}
	delete [] dblRES;

}
//---------------------------------------------------------------------------
int SolverLinear_GMRES::iluk(int *n, double *a, int *ja, int *ia, int *lfil,
                        double*& aluold, int*& jluold, int *ju, int*& levsold, int *iwk,
                        double *w, int *jw, int *ierr)
{
	/* System generated locals */
	int i__1, i__2, i__3, i__4;

	/* Local variables */
	static double fact;
	static int lenl, jlev, lenu, jpos, jrow, i__, j, k;
	static double s, t;
	static int j1, j2, n2, ii, jj, ju0;

	std::vector<double> alu((*iwk)+2,0.);
	std::vector<int> jlu((*iwk)+2,0);
	std::vector<int> levs((*iwk)+2,0);

/*=========================================================================*
 *     SPARSKIT ROUTINE ILUK -- ILU WITH LEVEL OF FILL-IN OF K (ILU(k))    *
 *=========================================================================*
 *                                                                         *
 * on entry:                                                               *
 * ==========                                                              *
 * n       = integer. The row dimension of the matrix A. The matrix        *
 * a,ja,ia = matrix stored in Compressed Sparse Row format.                *
 * lfil    = integer. The fill-in parameter. Each element whose            *
 *           leve-of-fill exceeds lfil during the ILU process is dropped.  *
 *           lfil must be .ge. 0                                           *
 * tol     = real*8. Sets the threshold for dropping small terms in the    *
 *           factorization. See below for details on dropping strategy.    *
 * iwk     = integer. The minimum length of arrays alu, jlu, and levs.     *
 *                                                                         *
 * On return:                                                              *
 * ===========                                                             *
 * alu,jlu = matrix stored in Modified Sparse Row (MSR) format containing  *
 *           the L and U factors together. The diagonal (stored in         *
 *           alu(1:n) ) is inverted. Each i-th row of the alu,jlu matrix   *
 *           contains the i-th row of L (excluding the diagonal entry=1)   *
 *           followed by the i-th row of U.                                *
 * ju      = integer array of length n containing the pointers to          *
 *           the beginning of each row of U in the matrix alu,jlu.         *
 * levs    = integer (work) array of size iwk -- which contains the        *
 *           levels of each element in alu, jlu.                           *
 * ierr    = integer. Error message with the following meaning.            *
 *           ierr  = 0    --> successful return.                           *
 *           ierr .gt. 0  --> zero pivot encountered at step number ierr.  *
 *           ierr  = -1   --> Error. input matrix may be wrong.            *
 *                            (The elimination process has generated a     *
 *                            row in L or U whose length is .gt.  n.)      *
 *           ierr  = -2   --> The matrix L overflows the array al.         *
 *           ierr  = -3   --> The matrix U overflows the array alu.        *
 *           ierr  = -4   --> Illegal value for lfil.                      *
 *           ierr  = -5   --> zero row encountered in A or U.              *
 *                                                                         *
 * work arrays:                                                            *
 * =============                                                           *
 * jw      = integer work array of length 3*n.                             *
 * w       = real work array of length n                                   *
 *                                                                         *
 * Notes/known bugs: This is not implemented efficiently storage-wise.     *
 *       For example: Only the part of the array levs(*) associated with   *
 *       the U-matrix is needed in the routine.. So some storage can       *
 *       be saved if needed. The levels of fills in the LU matrix are      *
 *       output for information only -- they are not needed by LU-solve.   *
 *                                                                         *
 *=========================================================================*
 *                                                                         *
 * w, ju (1:n) store the working array [1:ii-1 = L-part, ii:n = u]         *
 * jw(n+1:2n)  stores the nonzero indicator.                               *
 *                                                                         *
 * Notes:                                                                  *
 * ------                                                                  *
 * All the diagonal elements of the input matrix must be  nonzero.         *
 *                                                                         *
 *=========================================================================*/
/*     locals */
    /* Parameter adjustments */
	--jw;
	--w;
	--ju;
	--ia;
	--a;
	--ja;
//    --alu;
//    --jlu;
//    --levs;

   /* Function Body */
	if (*lfil < 0) {
		goto L998;
	}
/* ----------------------------------------------------------------------- */
/*     initialize ju0 (points to next element to be added to alu,jlu) */
/*     and pointer array. */
/* ----------------------------------------------------------------------- */
	n2 = *n + *n;
	ju0 = *n + 2;
	jlu[1] = ju0;

/*     initialize nonzero indicator array + levs array -- */

	i__1 = *n << 1;
	for (j = 1; j <= i__1; ++j) {
		jw[j] = 0;
/* L1: */
	}
/* ----------------------------------------------------------------------- */
/*     beginning of main loop. */
/* ----------------------------------------------------------------------- */
	i__1 = *n;
	for (ii = 1; ii <= i__1; ++ii) {
		j1 = ia[ii];
		j2 = ia[ii + 1] - 1;

/*     unpack L-part and U-part of row of A in arrays w */

		lenu = 1;
		lenl = 0;
		jw[ii] = ii;
		w[ii] = (float)0.;
		jw[*n + ii] = ii;

		i__2 = j2;
		for (j = j1; j <= i__2; ++j) {
			k = ja[j];
			t = a[j];
			if (t == (float)0.) {
				goto L170;
			}
			if (k < ii) {
				++lenl;
				jw[lenl] = k;
				w[lenl] = t;
				jw[n2 + lenl] = 0;
				jw[*n + k] = lenl;
			} 
			else if (k == ii) {
				w[ii] = t;
				jw[n2 + ii] = 0;
			} 
			else {
				++lenu;
				jpos = ii + lenu - 1;
				jw[jpos] = k;
				w[jpos] = t;
				jw[n2 + jpos] = 0;
				jw[*n + k] = jpos;
			}
		L170:
			;
		}

		jj = 0;

/*     eliminate previous rows */

		L150:
			++jj;
		if (jj > lenl) {
			goto L160;
		}
/* ----------------------------------------------------------------------- */
/*     in order to do the elimination in the correct order we must select */
/*     the smallest column index among jw(k), k=jj+1, ..., lenl. */
/* ----------------------------------------------------------------------- */
		jrow = jw[jj];
		k = jj;

/*     determine smallest column index */

		i__2 = lenl;
		for (j = jj + 1; j <= i__2; ++j) {
			if (jw[j] < jrow) {
				jrow = jw[j];
				k = j;
			}
			/* L151: */
		}

		if (k != jj) {
/*     exchange in jw */
			j = jw[jj];
			jw[jj] = jw[k];
			jw[k] = j;
/*     exchange in jw(n+  (pointers/ nonzero indicator). */
			jw[*n + jrow] = jj;
			jw[*n + j] = k;
/*     exchange in jw(n2+  (levels) */
			j = jw[n2 + jj];
			jw[n2 + jj] = jw[n2 + k];
			jw[n2 + k] = j;
/*     exchange in w */
			s = w[jj];
			w[jj] = w[k];
			w[k] = s;
		}

/*     zero out element in row by resetting jw(n+jrow) to zero. */

		jw[*n + jrow] = 0;

/*     get the multiplier for row to be eliminated (jrow) + its level */

		fact = w[jj] * alu[jrow];
		jlev = jw[n2 + jj];
		if (jlev > *lfil) {
			goto L150;
		}

/*     combine current row and row jrow */

		i__2 = jlu[jrow + 1] - 1;
		for (k = ju[jrow]; k <= i__2; ++k) {
			s = fact * alu[k];
			j = jlu[k];
			jpos = jw[*n + j];
			if (j >= ii) {

/*     dealing with upper part. */

				if (jpos == 0) {

/*     this is a fill-in element */

					++lenu;
					if (lenu > *n) {
						goto L995;
					}
					i__ = ii + lenu - 1;
					jw[i__] = j;
					jw[*n + j] = i__;
					w[i__] = -s;
					jw[n2 + i__] = jlev + levs[k] + 1;
				}
				else {

/*     this is not a fill-in element */

					w[jpos] -= s;
/* Computing MIN */
					i__3 = jw[n2 + jpos], i__4 = jlev + levs[k] + 1;
					jw[n2 + jpos] = std::min(i__3,i__4);
				}
			}
			else {

/*     dealing with lower part. */

				if (jpos == 0) {

/*     this is a fill-in element */

					++lenl;
					if (lenl > *n) {
						goto L995;
					}
					jw[lenl] = j;
					jw[*n + j] = lenl;
					w[lenl] = -s;
					jw[n2 + lenl] = jlev + levs[k] + 1;
				} 
				else {

/*     this is not a fill-in element */

					w[jpos] -= s;
/* Computing MIN */
					i__3 = jw[n2 + jpos], i__4 = jlev + levs[k] + 1;
					jw[n2 + jpos] = std::min(i__3,i__4);
				}
			}
			/* L203: */
		}
		w[jj] = fact;
		jw[jj] = jrow;
		goto L150;
		L160:

/*     reset double-pointer to zero (U-part) */

		i__2 = lenu;
		for (k = 1; k <= i__2; ++k) {
			jw[*n + jw[ii + k - 1]] = 0;
		/* L308: */
		}

/*     update l-matrix */

		i__2 = lenl;
		for (k = 1; k <= i__2; ++k) {

			if (ju0 > *iwk) {    
			
			//cout << "\tju0 = " << ju0 << "\tiwk = " << *iwk << "\tsize van alu = " << alu.size() << endl;
			
				*ierr = -2;
				//return 0;
				//goto L996;
				alu.push_back(0.);
				jlu.push_back(0);
				levs.push_back(0);
				*iwk += 1;
			}
			if (jw[n2 + k] <= *lfil) {
				alu[ju0] = w[k];
				jlu[ju0] = jw[k];
				++ju0;
			}
			/* L204: */
		}
	
	//cout << "einde van for lus : ju0 = " << ju0 << endl;

	//return 0;
	//goto L996;

/*     save pointer to beginning of row ii of U */

		ju[ii] = ju0;

/*     update u-matrix */

		i__2 = ii + lenu - 1;
		for (k = ii + 1; k <= i__2; ++k) {
			if (jw[n2 + k] <= *lfil) {
				if (ju0 > *iwk) {    
					alu.push_back(0.);
					jlu.push_back(0);
					levs.push_back(0);
					*iwk += 1;
				}

				jlu[ju0] = jw[k];
				alu[ju0] = w[k];
				levs[ju0] = jw[n2 + k];
				++ju0;
			}
		/* L302: */
		}
		if (w[ii] == (float)0.) {
			goto L999;
		}

		alu[ii] = 1. / w[ii];

/*     update pointer to beginning of next row of U. */

		jlu[ii + 1] = ju0;
/* ----------------------------------------------------------------------- */
/*     end main loop */
/* ----------------------------------------------------------------------- */
/* L500: */
    }

	delete [] aluold;
	delete [] jluold;
	//delete [] levsold;

	aluold = new double[*iwk+1]; 
	if (aluold == 0)
    {
           //cout << "Not enough memory to allocate lu4" << endl;
    }
	jluold = new int[*iwk+1]; 
	if (jluold == 0)
    {
           //cout << "Not enough memory to allocate lu5" << endl;
    }
	//levsold = new int[*iwk+1];

	for (unsigned i=0;i<*iwk;++i)
	{
		aluold[i] = alu[i+1];
		jluold[i] = jlu[i+1];
		//levsold[i] = levs[i+1];
	}

    *ierr = 0;
	return *iwk;
    //return 0;

/*     incomprehensible error. Matrix must be wrong. */

L995:
	*ierr = -1;
	return 0;

/*     insufficient storage in L. */

L996:
	*ierr = -2;
	return 0;

/*     insufficient storage in U. */

/* L997: */
	*ierr = -3;
	return 0;

// illegal lfil entered.
L998:
	*ierr = -4;
	return 0;

// zero row encountered in A or U.
L999:
	*ierr = -5;
	return 0;
}
//---------------------------------------------------------------------------
double SolverLinear_GMRES::ddot(int *n, double *dx, int *incx, double *dy, int *incy)
{
	/* System generated locals */
	int i__1;
	double ret_val;

	/* Local variables */
	static int i__, m;
	static double dtemp;
	static int ix, iy, mp1;


/*     forms the dot product of two vectors. */
/*     uses unrolled loops for increments equal to one. */
/*     jack dongarra, linpack, 3/11/78. */


    /* Parameter adjustments */
    --dy;
    --dx;

    /* Function Body */
    ret_val = 0.;
    dtemp = 0.;
    if (*n <= 0) {
	return ret_val;
    }
    if (*incx == 1 && *incy == 1) {
	goto L20;
    }

/*        code for unequal increments or equal increments */
/*          not equal to 1 */

    ix = 1;
    iy = 1;
    if (*incx < 0) {
	ix = (-(*n) + 1) * *incx + 1;
    }
    if (*incy < 0) {
	iy = (-(*n) + 1) * *incy + 1;
    }
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	dtemp += dx[ix] * dy[iy];
	ix += *incx;
	iy += *incy;
/* L10: */
    }
    ret_val = dtemp;
    return ret_val;

/*        code for both increments equal to 1 */


/*        clean-up loop */

L20:
    m = *n % 5;
    if (m == 0) {
	goto L40;
    }
    i__1 = m;
    for (i__ = 1; i__ <= i__1; ++i__) {
	dtemp += dx[i__] * dy[i__];
/* L30: */
    }
    if (*n < 5) {
	goto L60;
    }
L40:
    mp1 = m + 1;
    i__1 = *n;
    for (i__ = mp1; i__ <= i__1; i__ += 5) {
	dtemp = dtemp + dx[i__] * dy[i__] + dx[i__ + 1] * dy[i__ + 1] + dx[
		i__ + 2] * dy[i__ + 2] + dx[i__ + 3] * dy[i__ + 3] + dx[i__ + 
		4] * dy[i__ + 4];
/* L50: */
    }
L60:
    ret_val = dtemp;
    return ret_val;
}
//---------------------------------------------------------------------------
int SolverLinear_GMRES::daxpy(int *n, double *da, double *dx, int *incx, double *dy, int *incy)
{
    /* System generated locals */
    int i__1;

    /* Local variables */
    static int i__, m, ix, iy, mp1;


/*     constant times a vector plus a vector. */
/*     uses unrolled loops for increments equal to one. */
/*     jack dongarra, linpack, 3/11/78. */


    /* Parameter adjustments */
    --dy;
    --dx;

    /* Function Body */
    if (*n <= 0) {
	return 0;
    }
    if (*da == 0.) {
	return 0;
    }
    if (*incx == 1 && *incy == 1) {
	goto L20;
    }

/*        code for unequal increments or equal increments */
/*          not equal to 1 */

    ix = 1;
    iy = 1;
    if (*incx < 0) {
	ix = (-(*n) + 1) * *incx + 1;
    }
    if (*incy < 0) {
	iy = (-(*n) + 1) * *incy + 1;
    }
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	dy[iy] += *da * dx[ix];
	ix += *incx;
	iy += *incy;
/* L10: */
    }
    return 0;

/*        code for both increments equal to 1 */


/*        clean-up loop */

L20:
    m = *n % 4;
    if (m == 0) {
	goto L40;
    }
    i__1 = m;
    for (i__ = 1; i__ <= i__1; ++i__) {
	dy[i__] += *da * dx[i__];
/* L30: */
    }
    if (*n < 4) {
	return 0;
    }
L40:
    mp1 = m + 1;
    i__1 = *n;
    for (i__ = mp1; i__ <= i__1; i__ += 4) {
	dy[i__] += *da * dx[i__];
	dy[i__ + 1] += *da * dx[i__ + 1];
	dy[i__ + 2] += *da * dx[i__ + 2];
	dy[i__ + 3] += *da * dx[i__ + 3];
/* L50: */
    }
    return 0;
}
//---------------------------------------------------------------------------
double SolverLinear_GMRES::dnrm2(int *n, double *dx, int *incx)
{
    /* initialized datFile */

    static double zero = 0.;
    static double one = 1.;
    static double cutlo = 8.232e-11;
    static double cuthi = 1.304e19;

    /* Format strings */
    static char fmt_30[] = "";
    static char fmt_50[] = "";
    static char fmt_70[] = "";
    static char fmt_110[] = "";

    /* System generated locals */
    int i__1, i__2;
    double ret_val, d__1;

    /* Builtin functions */
    //double sqrt();

    /* Local variables */
    static double xmax;
    static int next, i__, j, nn;
    static double hitest, sum;

    /* Assigned format variables */
    static char *next_fmt;

    /* Parameter adjustments */
    --dx;

    /* Function Body */

/*     euclidean norm of the n-vector stored in dx() with storage */
/*     increment incx . */
/*     if    n .le. 0 return with result = 0. */
/*     if n .ge. 1 then incx must be .ge. 1 */


    if (*n > 0) {
	goto L10;
    }
    ret_val = zero;
    goto L300;

L10:
    next = 0;
    next_fmt = fmt_30;
    sum = zero;
    nn = *n * *incx;
/*                                                 begin main loop */
    i__ = 1;
L20:
    switch ((int)next) {
	case 0: goto L30;
	case 1: goto L50;
	case 2: goto L70;
	case 3: goto L110;
    }
L30:
    if ((d__1 = dx[i__], fabs(d__1)) > cutlo) {
	goto L85;
    }
    next = 1;
    next_fmt = fmt_50;
    xmax = zero;

/*                        phase 1.  sum is zero */

L50:
    if (dx[i__] == zero) {
	goto L200;
    }
    if ((d__1 = dx[i__], fabs(d__1)) > cutlo) {
	goto L85;
    }

/*                                prepare for phase 2. */
    next = 2;
    next_fmt = fmt_70;
    goto L105;

/*                                prepare for phase 4. */

L100:
    i__ = j;
    next = 3;
    next_fmt = fmt_110;
    sum = sum / dx[i__] / dx[i__];
L105:
    xmax = (d__1 = dx[i__], fabs(d__1));
    goto L115;

/*                   phase 2.  sum is small. */
/*                             scale to avoid destructive underflow. */

L70:
    if ((d__1 = dx[i__], fabs(d__1)) > cutlo) {
	goto L75;
    }

/*                     common code for phases 2 and 4. */
/*                     in phase 4 sum is large.  scale to avoid overflow. */

L110:
    if ((d__1 = dx[i__], fabs(d__1)) <= xmax) {
	goto L115;
    }
/* Computing 2nd power */
    d__1 = xmax / dx[i__];
    sum = one + sum * (d__1 * d__1);
    xmax = (d__1 = dx[i__], fabs(d__1));
    goto L200;

L115:
/* Computing 2nd power */
    d__1 = dx[i__] / xmax;
    sum += d__1 * d__1;
    goto L200;


/*                  prepare for phase 3. */

L75:
    sum = sum * xmax * xmax;


/*     for real or d.p. set hitest = cuthi/n */
/*     for complex      set hitest = cuthi/(2*n) */

L85:
    hitest = cuthi / (double) (*n);

/*                   phase 3.  sum is mid-range.  no scaling. */

    i__1 = nn;
    i__2 = *incx;
    for (j = i__; i__2 < 0 ? j >= i__1 : j <= i__1; j += i__2) {
	if ((d__1 = dx[j], fabs(d__1)) >= hitest) {
	    goto L100;
	}
/* L95: */
/* Computing 2nd power */
	d__1 = dx[j];
	sum += d__1 * d__1;
    }
    ret_val = sqrt(sum);
    goto L300;

L200:
    i__ += *incx;
    if (i__ <= nn) {
	goto L20;
    }

/*              end of main loop. */

/*              compute square root and adjust for scaling. */

    ret_val = xmax * sqrt(sum);
L300:
    return ret_val;
}
//---------------------------------------------------------------------------
void SolverLinear_GMRES::amux(int *n, double *x, double *y, double *a, int *ja, int *ia)
{
   /* System generated locals */
   int i__1, i__2;

   /* Local variables */
   static int i__, k;
   static double t;

/*=========================================================================*
 *         A times a vector                                                *
 *=========================================================================*
 * multiplies a matrix by a vector using the dot product form              *
 * Matrix A is stored in compressed sparse row storage.                    *
 *                                                                         *
 * on entry:                                                               *
 * ----------                                                              *
 * n     = row dimension of A                                              *
 * x     = real array of length equal to the column dimension of           *
 *         the A matrix.                                                   *
 * a, ja,                                                                  *
 *    ia = input matrix in compressed sparse row format.                   *
 *                                                                         *
 * on return:                                                              *
 * -----------                                                             *
 * y     = real array of length n, containing the product y=Ax             *
 *                                                                         *
 *=========================================================================*/
/* local variables */

/* ----------------------------------------------------------------------- */
   /* Parameter adjustments */
   --ia;
   --ja;
   --a;
   --y;
   --x;

   /* Function Body */
   i__1 = *n;
   for (i__ = 1; i__ <= i__1; ++i__) {

      /* compute the inner product of row i with vector x */
      t = 0.;
      i__2 = ia[i__ + 1] - 1;
      for (k = ia[i__]; k <= i__2; ++k) {
         t += a[k] * x[ja[k]];
      }

      /* store result in y(i) */
      y[i__] = t;
      
   }
}
//---------------------------------------------------------------------------
void SolverLinear_GMRES::lusol(int *n, double *y, double *x, double *alu, int *jlu, int *ju)
{
   /* System generated locals */
   int i__1, i__2;

   /* Local variables */
   static int i__, k;

/*=========================================================================*
 *                                                                         *
 * This routine solves the system (LU) x = y,                              *
 * given an LU decomposition of a matrix stored in (alu, jlu, ju)          *
 * modified sparse row format                                              *
 *                                                                         *
 *=========================================================================*
 * on entry:                                                               *
 * n   = dimension of system                                               *
 * y   = the right-hand-side vector                                        *
 * alu, jlu, ju                                                            *
 *     = the LU matrix as provided from the ILU routines.                  *
 *                                                                         *
 * on return                                                               *
 * x   = solution of LU x = y.                                             *
 *=========================================================================*
 *                                                                         *
 * Note: routine is in place: call lusol (n, x, x, alu, jlu, ju)           *
 *       will solve the system with rhs x and overwrite the result on x .  *
 *                                                                         *
 *=========================================================================*/
/* local variables */


/* forward solve */

   /* Parameter adjustments */
   --x;
   --y;
   --alu;
   --jlu;
   --ju;
   /* Function Body */
   i__1 = *n;
   for (i__ = 1; i__ <= i__1; ++i__) {
      x[i__] = y[i__];
      i__2 = ju[i__] - 1;
      for (k = jlu[i__]; k <= i__2; ++k) {
         x[i__] -= alu[k] * x[jlu[k]];
      }
   }

   /* backward solve. */
   for (i__ = *n; i__ >= 1; --i__) {
      i__1 = jlu[i__ + 1] - 1;
      for (k = ju[i__]; k <= i__1; ++k) {
         x[i__] -= alu[k] * x[jlu[k]];
      }
      x[i__] = alu[i__] * x[i__];
   }
}
//---------------------------------------------------------------------------
void SolverLinear_GMRES::pgmres(int *n, int *im, double *rhs, double *sol, double *vv,
             double *eps, int *maxits, int*iout,
             double *aa, int *ja, int *ia,
             double *alu, int *jlu, int *ju,
             int *ierr)
{
   /* initialized datFile */
   static double epsmac = 1e-16;

   /* System generated locals */
   int vv_dim1, vv_offset, i__1, i__2;
   double d__1, d__2;

   /* Local variables */
   static double c__[50];
   static int i__, j, k;
   static double s[50], t;
   static int i1, k1;
   static int n1;
   static double hh[2550]	/* was [51][50] */;
   static int ii, jj;
   static double ro, rs[51], gam;
   static int its;
   static double eps1;

   //static AnsiString Message;


/*=========================================================================*
 *                                                                         *
 *                 *** ILUT - Preconditioned GMRES ***                     *
 *                                                                         *
 *=========================================================================*
 * This is a simple version of the ILUT preconditioned GMRES algorithm.    *
 * The ILUT preconditioner uses a dual strategy for dropping elements      *
 * instead  of the usual level of-fill-in approach. See details in ILUT    *
 * subroutine documentation. PGMRES uses the L and U matrices generated    *
 * from the subroutine ILUT to precondition the GMRES algorithm.           *
 * The preconditioning is applied to the right. The stopping criterion     *
 * utilized is based simply on reducing the residual norm by epsilon.      *
 * This preconditioning is more reliable than ilu0 but requires more       *
 * storage. It seems to be much less prone to difficulties related to      *
 * strong nonsymmetries in the matrix. We recommend using a nonzero tol    *
 * (tol=.005 or .001 usually give good results) in ILUT. Use a large       *
 * lfil whenever possible (e.g. lfil = 5 to 10). The higher lfil the       *
 * more reliable the code is. Efficiency may also be much improved.        *
 * Note that lfil=n and tol=0.0 in ILUT  will yield the same factors as    *
 * Gaussian elimination without pivoting.                                  *
 *                                                                         *
 * ILU(0) and MILU(0) are also provided for comparison purposes            *
 * USAGE: first call ILUT or ILU0 or MILU0 to set up preconditioner and    *
 * then call pgmres.                                                       *
 *=========================================================================*
 * Coded by Y. Saad - This version dated May, 7, 1990.                     *
 *=========================================================================*
 * parameters                                                              *
 * -----------                                                             *
 * on entry:                                                               *
 * ==========                                                              *
 *                                                                         *
 * n     == integer. The dimension of the matrix.                          *
 * im    == size of krylov subspace:  should not exceed 50 in this         *
 *          version (can be reset by changing parameter command for        *
 *          kmax below)                                                    *
 * rhs   == real vector of length n containing the right hand side.        *
 *          Destroyed on return.                                           *
 * sol   == real vector of length n containing an initial guess to the     *
 *          solution on input. approximate solution on output              *
 * eps   == tolerance for stopping criterion. process is stopped           *
 *          as soon as ( ||.|| is the euclidean norm):                     *
 *          || current residual||/||initial residual|| <= eps              *
 * maxits== maximum number of iterations allowed                           *
 * iout  == output unit number number for printing intermediate results    *
 *          if (iout .le. 0) nothing is printed out.                       *
 *                                                                         *
 * aa, ja,                                                                 *
 * ia    == the input matrix in compressed sparse row format:              *
 *          aa(1:nnz)  = nonzero elements of A stored row-wise in order    *
 *          ja(1:nnz) = corresponding column indices.                      *
 *          ia(1:n+1) = pointer to beginning of each row in aa and ja.     *
 *          here nnz = number of nonzero elements in A = ia(n+1)-ia(1)     *
 *                                                                         *
 * alu,jlu== A matrix stored in Modified Sparse Row format containing      *
 *           the L and U factors, as computed by subroutine ilut.          *
 *                                                                         *
 * ju     == integer array of length n containing the pointers to          *
 *           the beginning of each row of U in alu, jlu as computed        *
 *           by subroutine ILUT.                                           *
 *                                                                         *
 * on return:                                                              *
 * ==========                                                              *
 * sol   == contains an approximate solution (upon successful return).     *
 * ierr  == integer. Error message with the following meaning.             *
 *          ierr = 0 --> successful return.                                *
 *          ierr = 1 --> convergence not achieved in itmax iterations.     *
 *          ierr =-1 --> the initial guess seems to be the exact           *
 *                       solution (initial residual computed was zero)     *
 *                                                                         *
 *=========================================================================*
 *                                                                         *
 * work arrays:                                                            *
 * =============                                                           *
 * vv    == work array of length  n x (im+1) (used to store the Arnoli     *
 *          basis)                                                         *
 *=========================================================================*
 * subroutines called :                                                    *
 * amux   : SPARSKIT routine to do the matrix by vector multiplication     *
 *          delivers y=Ax, given x  -- see SPARSKIT/BLASSM/amux            *
 * lusol : combined forward and backward solves (Preconditioning ope.)     *
 * BLAS1  routines.                                                        *
 *=========================================================================*
 *                                                                         *
 * arnoldi size should not exceed kmax=50 in this version..                *
 * to reset modify paramter kmax accordingly.                              *
 *=========================================================================*/
   /* Parameter adjustments */
   --ju;
   --ia;
   vv_dim1 = *n;
   vv_offset = 1 + vv_dim1 * 1;
   vv -= vv_offset;
   --sol;
   --rhs;
   --aa;
   --ja;
   --alu;
   --jlu;

   /* Function Body */
   n1 = *n + 1;
   its = 0;

   /* compute initial residual vector */
   amux(n, &sol[1], &vv[vv_offset], &aa[1], &ja[1], &ia[1]);
   i__1 = *n;
   for (j = 1; j <= i__1; ++j) {
      vv[j + vv_dim1] = rhs[j] - vv[j + vv_dim1];
   }

   /* outer loop starts here.. */
L20:
   ro = dnrm2(n, &vv[vv_offset], &c__1);
   if (*iout > 0 && its == 0) {
   }
   if (ro == 0.) {
     *ierr = -1;
     return;
	//goto L999;
   }
   t = 1. / ro;
   i__1 = *n;
   for (j = 1; j <= i__1; ++j) {
      vv[j + vv_dim1] *= t;
   }
   if (its == 0) {
      eps1 = *eps * ro;
   }

   /* initialize 1-st term  of rhs of hessenberg system.. */
   rs[0] = ro;
   i__ = 0;

L4:
   ++i__;
   ++its;
   i1 = i__ + 1;
   lusol(n, &vv[i__ * vv_dim1 + 1], &rhs[1], &alu[1], &jlu[1], &ju[1]);
   amux(n, &rhs[1], &vv[i1 * vv_dim1 + 1], &aa[1], &ja[1], &ia[1]);

   /* modified gram - schmidt... */
   i__1 = i__;
   for (j = 1; j <= i__1; ++j) {
      t = ddot(n, &vv[j * vv_dim1 + 1], &c__1, &vv[i1 * vv_dim1 + 1], &c__1);
      hh[j + i__ * 51 - 52] = t;
      d__1 = -t;
      daxpy(n, &d__1, &vv[j * vv_dim1 + 1], &c__1, &vv[i1 * vv_dim1 + 1], &c__1);
   }
   t = dnrm2(n, &vv[i1 * vv_dim1 + 1], &c__1);
   hh[i1 + i__ * 51 - 52] = t;
   if (t == 0.) {
      goto L58;
   }
   t = 1. / t;
   i__1 = *n;
   for (k = 1; k <= i__1; ++k) {
      vv[k + i1 * vv_dim1] *= t;
   }

   /* done with modified gram schimd and arnoldi step.. */
   /* now  update factorization of hh */
L58:
   if (i__ == 1) {
      goto L121;
   }

   /* perfrom previous transformations on i-th column of h */
   i__1 = i__;
   for (k = 2; k <= i__1; ++k) {
      k1 = k - 1;
      t = hh[k1 + i__ * 51 - 52];
      hh[k1 + i__ * 51 - 52] = c__[k1 - 1] * t + s[k1 - 1] * hh[k + i__ *51 - 52];
      hh[k + i__ * 51 - 52] = -s[k1 - 1] * t + c__[k1 - 1] * hh[k + i__ *51 - 52];
   }

L121:
   /* Computing 2nd power */
   d__1 = hh[i__ + i__ * 51 - 52];

   /* Computing 2nd power */
   d__2 = hh[i1 + i__ * 51 - 52];
   gam = sqrt(d__1 * d__1 + d__2 * d__2);

   /* if gamma is zero then any small value will do... */
   /* will affect only residual estimate */
   if (gam == 0.) {
      gam = epsmac;
   }

   /* get next plane rotation */
   c__[i__ - 1] = hh[i__ + i__ * 51 - 52] / gam;
   s[i__ - 1] = hh[i1 + i__ * 51 - 52] / gam;
   rs[i1 - 1] = -s[i__ - 1] * rs[i__ - 1];
   rs[i__ - 1] = c__[i__ - 1] * rs[i__ - 1];

   /* determine residual norm and test for convergence- */
   hh[i__ + i__ * 51 - 52] = c__[i__ - 1] * hh[i__ + i__ * 51 - 52]
                             + s[i__ - 1] * hh[i1 + i__ * 51 - 52];
   ro = (d__1 = rs[i1 - 1], fabs(d__1));

   /*if (*iout > 0) {
      Message = "Iteration ";
      Message += its;

      Message += "\t residual =  ";

      Message += ro;

      SolverForm->SolverOutput->Lines->Append(Message);

   }*/

   //unsigned remainder = (100u*its)%(*maxits);
   //if (remainder == 0)
   //  SolverForm->ProgressBar1->Position = 100u*its/(*maxits);

   if (i__ < *im && ro > eps1) {
      goto L4;
   }

   /* now compute solution. first solve upper triangular system. */
   rs[i__ - 1] /= hh[i__ + i__ * 51 - 52];
   i__1 = i__;
   for (ii = 2; ii <= i__1; ++ii) {
      k = i__ - ii + 1;
      k1 = k + 1;
      t = rs[k - 1];
      i__2 = i__;
      for (j = k1; j <= i__2; ++j) {
	 t -= hh[k + j * 51 - 52] * rs[j - 1];
      }
      rs[k - 1] = t / hh[k + k * 51 - 52];
   }

   /* form linear combination of v(*,i)'s to get solution */
   t = rs[0];
   i__1 = *n;
   for (k = 1; k <= i__1; ++k) {
      rhs[k] = vv[k + vv_dim1] * t;
   }
   i__1 = i__;
   for (j = 2; j <= i__1; ++j) {
      t = rs[j - 1];
      i__2 = *n;
      for (k = 1; k <= i__2; ++k) {
	 rhs[k] += t * vv[k + j * vv_dim1];
      }
   }

   /* call preconditioner. */
   lusol(n, &rhs[1], &rhs[1], &alu[1], &jlu[1], &ju[1]);
   i__1 = *n;
   for (k = 1; k <= i__1; ++k) {
      sol[k] += rhs[k];
   }

   /* restart outer loop  when necessary */
   if (ro <= eps1) {
      goto L990;
   }
   if (its >= *maxits) {
      goto L991;
   }

   /* else compute residual vector and continue.. */
   i__1 = i__;
   for (j = 1; j <= i__1; ++j) {
      jj = i1 - j + 1;
      rs[jj - 2] = -s[jj - 2] * rs[jj - 1];
      rs[jj - 1] = c__[jj - 2] * rs[jj - 1];
   }
   i__1 = i1;
   for (j = 1; j <= i__1; ++j) {
      t = rs[j - 1];
      if (j == 1) {
	 t += -1.;
      }
      daxpy(n, &t, &vv[j * vv_dim1 + 1], &c__1, &vv[vv_offset], &c__1);
   }

   /* restart outer loop. */
    goto L20;
    
L990:
    *ierr = 0;

L991:
    *ierr = 1;
    
L999:
    *ierr = -1;
}
//---------------------------------------------------------------------------
void SolverLinear_GMRES::print(const std::string &filename)
{
  std::ofstream CHECK;
	CHECK.open(filename.c_str(), std::ofstream::out | std::ofstream::app);
	CHECK.precision(12);
	for (unsigned r = 0; r < nRows; r++)
	{
		unsigned iNonzero = 0;
		unsigned nNonzeros = Rows[r + 1] - Rows[r];
		for (unsigned c = 0; c < nRows; c++)
		{		
			if (iNonzero < nNonzeros)
			{
				unsigned c1 = Columns[Rows[r] - 1 + iNonzero] - 1;
				if (c == c1)
				{
					CHECK << A[getIndex(r, c)] << '\t';
					iNonzero++;
				}
				else 
				{
					CHECK << "~\t";
				}
			}			
			else 
			{
				CHECK << "~\t";
			}
		}
		CHECK << '\n';
	}
	CHECK.close();
	
	
	
	
	/*std::ofstream CHECK;
	CHECK.open(filename.c_str(), std::ofstream::out | std::ofstream::app);
	CHECK.precision(12);
	CHECK << "IA\n";
	for (unsigned r=0; r<nRows+1; r++)
	{
		CHECK << r << '\t' << Rows[r] << std::endl;
	}
	CHECK << "JA\n";
	for (unsigned c=0; c<N*nVariables*nVariables; c++)
	{
		CHECK << c << '\t' << Columns[c] << std::endl;
	}
	CHECK << "A\n";
	for (unsigned a=0; a<N*nVariables*nVariables; a++)
	{
		CHECK << a << '\t' << A[a] << std::endl;
	}
	CHECK << std::endl;
	CHECK.close();*/
	
	
	/*output << "Rows :" << endl;
   for (unsigned r=0; r<nNodes*(nIons+1)+1; r++) {
      output << Rows[r] << "\t";
   }
   output << "\n\nColumns :" << endl;
   for (unsigned c=0; c<N*(nIons+1)*(nIons+1); c++) {
      output << Columns[c] << "\t";
   }
   output << "\n\nA :" << endl;
   for (unsigned c=0; c<N*(nIons+1)*(nIons+1); c++) {
      output << A[c] << "\t";
   }
   output << "\n" << endl;
   unsigned index = 0;
   for (unsigned m=0; m<nNodes; m++) {
      for (unsigned i=0; i<nIons+1; i++) {
         for (unsigned n=0; n<nNodes; n++) {
            for (unsigned j=0; j<nIons+1; j++) {
               if (index == var(m,i)*nNodes*(nIons+1)+Columns[GetIndex(var(m,i),var(n,j))]-1)
                  output << A[GetIndex(var(m,i),var(n,j))] << "\t";
               else
                  output << "0" << "\t";
               index++;
            }
         }
         output << endl;
      }
   }
   output << endl;*/
}
//---------------------------------------------------------------------------