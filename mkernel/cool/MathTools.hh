#ifndef cool_MathTools_hh
#define cool_MathTools_hh

#include <cmath>
#include <vector>


//////////////////////////////////////////////////////////////////////////////


namespace MathTools {
  namespace MathFunctions {


double innerProd(const std::vector< double >& v1, const std::vector< double >& v2, const double& v=double());
void crossProd(const std::vector< double >& v1, const std::vector< double >& v2, std::vector< double >& r);


  }
}


//////////////////////////////////////////////////////////////////////////////


class RealVector : public std::vector< double >
{
 public:
  RealVector(const RealVector& other, bool dummy=false) : std::vector< double >(other) {}
  RealVector(const double& _v, const unsigned& _size) :
    std::vector< double >(_size,_v) {}

  // unary minus operator
  const RealVector& operator-() {
    for (std::size_t i=0; i<size(); ++i)
      operator[](i) = - operator[](i);
    return *this;
  }

  // binary operator on RealVector implementations, returning new RealVector
#define REALVECTOR_OP_REALVECTOR(__op__) \
  RealVector operator __op__ (const RealVector& other) { \
    const std::size_t s = std::min(size(),other.size()); \
    RealVector r(0.,s); \
    for (std::size_t i=0; i<s; ++i) \
      r[i] = operator[](i) __op__ other[i]; \
    return r; \
  }
REALVECTOR_OP_REALVECTOR(-)
REALVECTOR_OP_REALVECTOR(+)
REALVECTOR_OP_REALVECTOR(*)
REALVECTOR_OP_REALVECTOR(/)
#undef REALVECTOR_OP_REALVECTOR

  // binary operator on RealVector implementations, returning this RealVector
#define REALVECTOR_OP_REALVECTOR(__op__) \
  const RealVector& operator __op__ (const RealVector& v) { \
    const std::size_t s = std::min(size(),v.size()); \
    for (std::size_t i=0; i<s; ++i) \
      operator[](i) __op__ v[i]; \
    return *this; \
  }
REALVECTOR_OP_REALVECTOR(=)
REALVECTOR_OP_REALVECTOR(-=)
REALVECTOR_OP_REALVECTOR(+=)
REALVECTOR_OP_REALVECTOR(*=)
REALVECTOR_OP_REALVECTOR(/=)
#undef REALVECTOR_OP_REALVECTOR

  // binary operator on double implementations, returning new RealVector
#define REALVECTOR_OP_REAL(__op__) \
  RealVector operator __op__ (const double& v) { \
    RealVector r(0.,size()); \
    for (std::size_t i=0; i<size(); ++i) \
      r[i] = operator[](i) __op__ v; \
    return r; \
  }
REALVECTOR_OP_REAL(-)
REALVECTOR_OP_REAL(+)
REALVECTOR_OP_REAL(*)
REALVECTOR_OP_REAL(/)
#undef REALVECTOR_OP_REAL

  // binary operator on double implementations, returning this RealVector
#define REALVECTOR_OP_REAL(__op__) \
  const RealVector& operator __op__ (const double& v) { \
    for (std::size_t i=0; i<size(); ++i) \
      operator[](i) __op__ v; \
    return *this; \
  }
REALVECTOR_OP_REAL(=)
REALVECTOR_OP_REAL(-=)
REALVECTOR_OP_REAL(+=)
REALVECTOR_OP_REAL(*=)
REALVECTOR_OP_REAL(/=)
#undef REALVECTOR_OP_REAL

  double sqrNorm() { return MathTools::MathFunctions::innerProd(*this,*this); }

};


//////////////////////////////////////////////////////////////////////////////


#endif

