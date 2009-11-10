#ifndef mpoint_h
#define mpoint_h


#include <iostream>
#include <vector>
#include <cmath>
#include <cstdlib>


namespace m {


// forward declarations
struct mpoint;
#define OP(_O_)\
mpoint operator _O_ (const mpoint& a, const mpoint& b);
OP(+)
OP(-)
OP(*)
OP(/)
#undef OP
#define OP(_O_)\
mpoint operator _O_ (const mpoint& a, double d);
OP(+)
OP(-)
OP(*)
OP(/)
#undef OP
#define OP(_O_)\
mpoint operator _O_ (double d, const mpoint& a);
OP(+)
OP(-)
OP(*)
OP(/)
#undef OP
mpoint operator^(const mpoint& a, const mpoint& b);
mpoint operator-(const mpoint& a);


// cartesian point/vector
struct mpoint {
  double x,y,z;
  mpoint() { x=0.;  y=0.;  z=0.; }
  mpoint(double _x, double _y, double _z) { x=_x;  y=_y;  z=_z; }
  //mpoint(const mpoint& o) { x=o.x;  y=o.y;  z=o.z; }
#define OP(_O_) \
  mpoint& operator _O_ (const mpoint& o) { x _O_ o.x;  y _O_ o.y;  z _O_ o.z;  return *this; }
  OP(+=)
  OP(-=)
  OP(*=)
  OP(/=)
#undef OP
#define OP(_O_)\
  mpoint& operator _O_ (double d) { x _O_ d;  y _O_ d;  z _O_ d;  return *this; }
  OP(+=)
  OP(-=)
  OP(*=)
  OP(/=)
#undef OP

  // internal and external products
  double operator&(const mpoint& o) const {
    mpoint t(*this);
    t *= o;
    return t.x+t.y+t.z;
  }
  mpoint& operator^=(const mpoint& b) {
    mpoint a(*this);
    x=a.y*b.z-a.z*b.y;  y=a.z*b.x-a.x*b.z;  z=a.x*b.y-a.y*b.x;
    return *this;
  }

  // usual norms, to use as a point
  double length(const mpoint& o) const { return sqrt(norm2(o)); }
  double norm2(const mpoint& o) const  { return (x-o.x)*(x-o.x) + (y-o.y)*(y-o.y) + (z-o.z)*(z-o.z); }

  // usual norms, to use as a vector
  double length() const { return sqrt(norm2()); }
  double norm2() const  { return x*x+y*y+z*z; }
  mpoint& normalize() { return operator/=(length()); }

  // geometric utilities
  bool intriangle(const mpoint& p1, const mpoint& p2, const mpoint& p3) const
  {
    // NOTE: this is actually a 2d check, assuming z=0.
    // vectors, dot products, barycentric coordinates and check
    const mpoint v0(p1-p3);
    const mpoint v1(p2-p3);
    const mpoint v2((*this)-p3);
    const double dot00 = v0&v0;
    const double dot01 = v0&v1;
    const double dot02 = v0&v2;
    const double dot11 = v1&v1;
    const double dot12 = v1&v2;
    const double f = 1./(dot00*dot11 - dot01*dot01);
    const double u = (dot11*dot02 - dot01*dot12) * f;
    const double v = (dot00*dot12 - dot01*dot02) * f;
    return (u>=0) && (v>=0) && (u+v<=1);
  }

};


// cartesian point/vector operators
std::ostream& operator<<(std::ostream& o, const mpoint& p) {
#define SIMPLE(_X_) (std::abs(_X_)>1.e-9? _X_ : 0.)
  o << ' ' << SIMPLE(p.x) << ' ' << SIMPLE(p.y) << ' ' << SIMPLE(p.z) << ' ';
#undef SIMPLE
  return o;
}
#define OP(_O_)\
mpoint operator _O_ (const mpoint& a, const mpoint& b) { mpoint r; r.x=a.x _O_ b.x; r.y=a.y _O_ b.y; r.z=a.z _O_ b.z; return r; }
OP(+)
OP(-)
OP(*)
OP(/)
#undef OP
#define OP(_O_)\
mpoint operator _O_ (const mpoint& a, double d) { mpoint r; r.x=a.x _O_ d; r.y=a.y _O_ d; r.z=a.z _O_ d; return r; }
OP(+)
OP(-)
OP(*)
OP(/)
#undef OP
#define OP(_O_)\
mpoint operator _O_ (double d, const mpoint& a) { mpoint r; r.x=a.x _O_ d; r.y=a.y _O_ d; r.z=a.z _O_ d; return r; }
OP(+)
OP(-)
OP(*)
OP(/)
#undef OP
mpoint operator^(const mpoint& a, const mpoint& b) { mpoint r(a); r^=b; return r; }
mpoint operator-(const mpoint& a) { mpoint r(a); r*=-1.; return r; }


// geometric utilities
inline int randomi(int min=0, int max=1)
{
  const double r = (double) max - (double) min + 1.;
  return min + (int)(r*rand()/(RAND_MAX+1.));
}
inline double randomd(double min=0., double max=1.)
{
  const int r = randomi(0,100000);
  return (double)(r/100000.) * (max-min) + min;
}
mpoint intriangle(const mpoint& p1, const mpoint& p2, const mpoint& p3)
{
  double f1 = 10.;
  double f2 = 10.;
  while (f1+f2>1.) {
    f1 = randomd(0.,1.);
    f2 = randomd(0.,1.);
  }
  return mpoint(p3 + (p1-p3)*f1 + (p2-p3)*f2);
}
mpoint inlinesegment(const mpoint& p1, const mpoint& p2)
{
  return p2 + (p1-p2)*randomd(0.,1.);
}
std::vector< double > angles(const mpoint& p0, const mpoint& p1, const mpoint& p2)
{
  // determine triangle internal ([0;pi]) angles at each node
  std::vector< double > r(3);
  mpoint v01 = mpoint(p1-p0).normalize();
  mpoint v02 = mpoint(p2-p0).normalize();
  mpoint v12 = mpoint(p2-p1).normalize();
  r[0] = acos(v01&v02);
  r[1] = acos(v12&(-v01));
  r[2] = acos((-v02)&(-v12));
  return r;
}


}  // namespace m


#endif

