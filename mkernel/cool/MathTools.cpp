
#include "cool/MathTools.hh"
using namespace std;


namespace MathTools {
  namespace MathFunctions {


double innerProd(const vector< double >& v1, const vector< double >& v2, const double& v)
{
  double r(v);
  for (unsigned i=0; i<v1.size() && i<v2.size(); ++i)
    r += v1[i]*v2[i];
  return r;
}


void crossProd(const vector< double >& v1, const vector< double >& v2, vector< double >& r)
{
  if (v1.size()!=3 || v2.size()!=3)
    throw 42;
  r.assign(3,0.);
  r[0] =  v1[1]*v2[2]-v2[1]*v1[2];
  r[1] = -v1[0]*v2[2]+v2[0]*v1[2];
  r[2] =  v1[0]*v2[1]-v2[0]*v1[1];
}


  }
}

