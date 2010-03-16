//---------------------------------------------------------------------------

#include "MigrationTerm_1D_ISG.h"

//---------------------------------------------------------------------------


//---------------------------------------------------------------------------
MigrationTerm_1D_ISG::MigrationTerm_1D_ISG(unsigned nDimensions_, unsigned nNodes_, unsigned nVariables_, MITReM* mitrem_, ElementProps* elementProps_) 
	: MigrationTerm(nDimensions_, nNodes_,nVariables_, mitrem_, elementProps_)
{
}
//---------------------------------------------------------------------------
MigrationTerm_1D_ISG::~MigrationTerm_1D_ISG()
{
}
//---------------------------------------------------------------------------
void MigrationTerm_1D_ISG::calcMat(EmptyDoubleMatrix elementMat, DoubleVectorList coordinates, DoubleVectorList velocities, DoubleVectorList concentrations, DoubleList potentials, DoubleList temperatures, DoubleList densities, DoubleList voidFractions, DoubleVectorList magneticFieldVectors)
{
	double voidFraction = 0.;
	for (unsigned m=0; m<nNodes; m++)
	{
		voidFraction += voidFractions[m];
	}
	voidFraction *= 0.5;
	double Bruggeman = pow(1.-voidFraction,1.5);
	
	// Calculate coefficients
	for (unsigned m=0; m<nNodes; m++) 
	{
		mitrem->init(concentrations[m],potentials[m],temperatures[m],densities[m]);
		normal = elementProps->calcNormal(m,coordinates);
		for (unsigned c=0; c<nDimensions; c++)
		{
			normals[m][c] = normal[c];
		}
/*		for (unsigned i=0; i<nIons; i++) 
		{
			W[m][i] = mitrem->calcTransportMigrationFactor(i)*Bruggeman;
		}*/
	}
	double A, B;
	A=0; B=0;
	A=1e-4; B=1e-10;
/*	for (unsigned i=0; i<nIons; i++)
	{
		unsigned zi = mitrem->getIonChargeNumber(i);
		double Di = mitrem->getIonDiffusionConstant(i);
		A += zi*zi*Di*(concentrations[1][i] - concentrations[0][i])/elementSize;
		B += zi*zi*Di*(concentrations[0][i]*coordinates[1][0] - concentrations[1][i]*coordinates[0][0])/elementSize;
	}*/
	elementSize = elementProps->calcSize(coordinates);
//	double elementSize6 = 6.*elementSize;
/*	
	// Add to element matrix
	for (unsigned i=0; i<nIons; i++) 
	{
		double wi = mitrem->calcTransportMigrationFactor(i)*Bruggeman;
		for (unsigned m=0; m<nNodes; m++) 
		{
			unsigned eqmi = eq(m,i);
			if (A>1e-5)
			{
				Mmi[0] = -normals[m][0]*wi*(-A*(concentrations[0][i] - concentrations[1][i])*(coordinates[0][0] - coordinates[1][0])
																	  +(B*concentrations[0][i] - B*concentrations[1][i] - A*concentrations[1][i]*coordinates[0][0] + A*concentrations[0][i]*coordinates[1][0])
																 		 *(log(A*coordinates[0][0] + B) - log(A*coordinates[1][0] + B))
																 	 )
																	 /(A*elementSize*elementSize*(log(A*coordinates[0][0] + B) - log(A*coordinates[1][0] + B)));
//			Mmi[0] = ((2.*W[0][i] + W[1][i])*concentrations[0][i] + (W[0][i] + 2.*W[1][i])*concentrations[1][i])*normals[m][0]/elementSize6;
			}
			else
			{
				Mmi[0] = normals[m][0]*wi*(concentrations[0][i] + concentrations[1][i])/(2*(coordinates[0][0]-coordinates[1][0]));
			}
			for (unsigned n=0; n<nNodes; n++) 
			{
				elementMat[eqmi][var(n,nIons)] += normals[n][0]*Mmi[0];
			}
		}*/

	for (unsigned i=0; i<nIons; i++) 
	{
		double wi = mitrem->calcTransportMigrationFactor(i)*Bruggeman;
		for (unsigned m=0; m<nNodes; m++) 
		{
			unsigned eqmi = eq(m,i);
			if (A>1e-5)
			{
				Mmi[0] = -normals[m][0]*wi*(-A*(concentrations[0][i] - concentrations[1][i])*(coordinates[0][0] - coordinates[1][0])
																	  +(B*concentrations[0][i] - B*concentrations[1][i] - A*concentrations[1][i]*coordinates[0][0] + A*concentrations[0][i]*coordinates[1][0])
																 		 *(log(A*coordinates[0][0] + B) - log(A*coordinates[1][0] + B))
																 	 )
																	 /(A*elementSize*elementSize*(log(A*coordinates[0][0] + B) - log(A*coordinates[1][0] + B)));
//			Mmi[0] = ((2.*W[0][i] + W[1][i])*concentrations[0][i] + (W[0][i] + 2.*W[1][i])*concentrations[1][i])*normals[m][0]/elementSize6;
			}
			else
			{
				Mmi[0] = normals[m][0]*wi*(concentrations[0][i] + concentrations[1][i])/(2*(coordinates[0][0]-coordinates[1][0]));
			}
			for (unsigned n=0; n<nNodes; n++) 
			{
				elementMat[eqmi][var(n,nIons)] += normals[n][0]*Mmi[0];
			}
		}
	}
}
//---------------------------------------------------------------------------
void MigrationTerm_1D_ISG::calcJac(EmptyDoubleMatrix elementJac, DoubleVectorList coordinates, DoubleVectorList velocities, DoubleVectorList concentrations, DoubleList potentials, DoubleList temperatures, DoubleList densities, DoubleList voidFractions, DoubleVectorList magneticFieldVectors)
{
	double voidFraction = 0.;
	for (unsigned m=0; m<nNodes; m++)
	{
		voidFraction += voidFractions[m];
	}
	voidFraction *= 0.5;
	double Bruggeman = pow(1.-voidFraction,1.5);
	
	// Calculate coefficients
	for (unsigned c=0; c<nDimensions; c++)
	{
		gradU[c] = 0;
	}
	for (unsigned m=0; m<nNodes; m++) 
	{
		normal = elementProps->calcNormal(m,coordinates);
/*		for (unsigned c=0; c<nDimensions; c++)
		{
			gradU[c] += normal[c]*potentials[m];
		}*/
		mitrem->init(concentrations[m],potentials[m],temperatures[m],densities[m]);
		for (unsigned i=0; i<nIons; i++) 
		{
			W[m][i] = mitrem->calcTransportMigrationFactor(i)*Bruggeman;
		}
	}
	double A, B;
	A=0; B=0;
	for (unsigned i=0; i<nIons; i++)
	{
		unsigned zi = mitrem->getIonChargeNumber(i);
		double Di = mitrem->getIonDiffusionConstant(i);
		A += zi*zi*Di*(concentrations[1][i] - concentrations[0][i])/elementSize;
		B += zi*zi*Di*(concentrations[0][i]*coordinates[1][0] - concentrations[1][i]*coordinates[0][0])/elementSize;
	}
	elementSize = elementProps->calcSize(coordinates);
	// double elementSize6 = 6.*elementSize;
	for (unsigned c=0; c<nDimensions; c++)
	{
		gradU[c] /= elementSize;
	}
	
	// Add to element jacobian
	for (unsigned i=0; i<nIons; i++) 
	{
		// double wi = mitrem->calcTransportMigrationFactor(i)*Bruggeman;
		for (unsigned m=0; m<nNodes; m++) 
		{
			// unsigned eqmi = eq(m,i);
			unsigned p = (m+1)%2;
			Wmi[m][0] = (2.*W[m][i] + W[p][i])*normals[m][0]/6.;
			Wmi[p][0] = (W[m][i] + 2.*W[p][i])*normals[m][0]/6.;
/*			for (unsigned n=0; n<nNodes; n++) 
			{
				//elementJac[eqmi][var(n,i)] -= gradU[0]*Wmi[n][0];
				if (A>1e-5)
				{
					elementJac[eqmi][var(n,i)] += pow(double(-1.),int(n))*normals[m][0]*(potentials[0] - potentials[1])*wi
																				*(A*(coordinates[1][0] - coordinates[0][0]) 
																					+ (A*coordinates[1][0] + B)*(log(A*coordinates[0][0] + B) - log(A*coordinates[1][0] + B))
																				 )
																				 /(A*elementSize*elementSize*(log(A*coordinates[0][0] + B) - log(A*coordinates[1][0] + B)));
				}
				else
				{
					elementJac[eqmi][var(n,i)] -= normals[m][0]*(potentials[0] - potentials[1])*wi
																				/(2*(coordinates[0][0] - coordinates[1][0]));
				}
			}*/
		}
	}
}
//---------------------------------------------------------------------------

