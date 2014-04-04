#ifndef t_geo_h
#define t_geo_h

#include "mkernel.h"

// geometry checking transform module
class t_geo : public m::mtransform {
 public:
  t_geo() : eps(1.e-16) {}
  void transform(GetPot& o, m::mmesh& m);
  // quick check simplex element correctness outside of this module
  static bool cellgeom_ok(const std::vector< std::vector< double > >& coords);
 private:
  bool ispointintriangle(const double x, const double y, const double x1, const double y1, const double x2, const double y2, const double x3, const double y3);
  bool ispointintetrahedron(const double x, const double y, const double z, const double x1, const double y1, const double z1, const double x2, const double y2, const double z2, const double x3, const double y3, const double z3, const double x4, const double y4, const double z4);
  double eps;
 private:  // (interface from COOLFluiD)
  // computes scaled inward normals and volume for an element, following face
  // index defined by opposite node index
  void cellgeom(
    m::mmesh& m,
    int ic,
    struct local_node_struct *No_loc,
    double *vol,
    int *inc_min );
  // 3d vector product
  void vecprd3(double *v1, double *v2, double *v);
};


// (interface from COOLFluiD)
struct local_node_struct {
  int node;        // global node number
  double W[10];    // solution vector at nodes
  double Res[10];  // residual vector at nodes
  double norm[3];  // scaled inwards normal
  double C[4];     // scalar coefficients for cell
  double norm2;    // square of normal
};

#endif


