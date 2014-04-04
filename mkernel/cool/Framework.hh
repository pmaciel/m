
#ifndef cool_Framework_hh
#define cool_Framework_hh

#include <string>
#include "cool/MathTools.hh"


//////////////////////////////////////////////////////////////////////////////


#define CFNULL ((void*) 0)
enum CFDim { DIM_0D=0, DIM_1D=1, DIM_2D=2, DIM_3D=3 };
enum CoordXYZ { XX=0, YY=1, ZZ=2 };

typedef unsigned CFuint;
typedef double CFreal;

void cf_always_assert_desc(const std::string& m, const bool a);
void cf_assert_desc(const std::string& m, const bool a);
void cf_assert(const bool a);


//////////////////////////////////////////////////////////////////////////////


namespace COOLFluiD {
  namespace Framework {


typedef RealVector Node;
typedef RealVector State;


  }
}


//////////////////////////////////////////////////////////////////////////////


#endif

