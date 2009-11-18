#ifndef t_math_h
#define t_math_h

#include "mkernel.h"

// mathematical transform module
class t_math : public m::mtransform {
 public:
  void transform(GetPot& o, m::mmesh& m);
 private:
  void rotate(std::vector< double >& vx, std::vector< double >& vy, std::vector< std::vector< double >* >& vvx, std::vector< std::vector< double >* >& vvy, const double a);
  void scale(std::vector< double >& v, const double f);
  void translate(std::vector< double >& v, const double d);
};

#endif


