
#include <sstream>
#include <string>
#include "samg.h"


namespace m {
  namespace aux {


struct s_parm {

  // construtors for hidden/secondary parameters
  s_parm(const std::string& _n, SAMG_C_CALLCONV (*_fpi)(int*),     const int&         _i) : n(_n), fpi(_fpi), fpd(NULL), fps(NULL), pd(NULL), pi(NULL), d(0.), i(_i), s(""), p1(0), p2(0) {}
  s_parm(const std::string& _n, SAMG_C_CALLCONV (*_fpd)(double*),  const double&      _d) : n(_n), fpi(NULL), fpd(_fpd), fps(NULL), pd(NULL), pi(NULL), d(_d), i(0),  s(""), p1(0), p2(0) {}
  s_parm(const std::string& _n, SAMG_C_CALLCONV (*_fps)(int*,int*),const std::string& _s) : n(_n), fpi(NULL), fpd(NULL), fps(_fps), pd(NULL), pi(NULL), d(0.), i(0),  s(_s), p1(0), p2(0) {}

  // constructors for primary double/integer parameters and defaults
  s_parm(const std::string& _n, double *_pd, const double& _d) : n(_n), fpi(NULL), fpd(NULL), fps(NULL), pd(_pd),  pi(NULL), d(_d), i(0),  s(""), p1(0), p2(0) {}
  s_parm(const std::string& _n, int    *_pi, const int&    _i) : n(_n), fpi(NULL), fpd(NULL), fps(NULL), pd(NULL), pi(_pi),  d(0.), i(_i), s(""), p1(0), p2(0) {}

  // constructor for integer sub-parameters
  s_parm(const std::string& _n, int    *_pi, int _p1, int _p2) : n(_n), fpi(NULL), fpd(NULL), fps(NULL), pd(NULL), pi(_pi),  d(0.), i(0),  s(""), p1(_p1), p2(_p2), radix(10) {}

  // parameter string conversion, setting/default and getting
  template< class T > T cvt(const std::string& s, T init=T()) const
  {
    T r(init);
    std::istringstream iss(s);
    iss >> r;
    return r;
  }
  void set(const std::string& s);
  void def();
  bool ismain() { return (p1? false:true); }
  std::string get();

  // parameter name
  std::string n;

  // pointers to functions for hidden/secondary parameters
  SAMG_C_CALLCONV (*fpi)(int*);
  SAMG_C_CALLCONV (*fpd)(double*);
  SAMG_C_CALLCONV (*fps)(int*,int*);

  // pointers to primary parameters and default values
  double *pd;
  int    *pi;
  double      d;
  int         i;
  std::string s;

  // integer parameter digit range begin/end and number radix
  int p1, p2, radix;

};


  }  // namespace aux
}  // namespace m

