
#include "mfactory.h"
#include "t_geo.h"

using namespace std;
using namespace m;


Register< mtransform,t_geo > mt_geo("-tg","geometry check (output to \"m2m_geometry.plt\")");


void t_geo::transform(GetPot& o, mmesh& m, const XMLNode& x)
{
  // number of positive/negative size and normals for reporting
  unsigned n_positivesize = 0;
  unsigned n_negativesize = 0;

  const unsigned dim = m.d();
  double vol  = 0.;
  int inc_min = 0;
  vector< local_node_struct > No_loc(dim+1);

  // each of these tetrahedra is one node swapped to check if splitting is correct
  /*
  swap(m.vz[0].e2n[10].n[0],m.vz[0].e2n[10].n[1]);
  swap(m.vz[0].e2n[11].n[0],m.vz[0].e2n[11].n[2]);
  swap(m.vz[0].e2n[12].n[0],m.vz[0].e2n[12].n[3]);
  swap(m.vz[0].e2n[13].n[1],m.vz[0].e2n[13].n[2]);
  swap(m.vz[0].e2n[14].n[1],m.vz[0].e2n[14].n[3]);
  swap(m.vz[0].e2n[15].n[2],m.vz[0].e2n[15].n[3]);
  */

  // output file
  ofstream f("m2m_geometry.plt",ios::trunc);
  f << "TITLE = \"problems\"" << endl
    << "VARIABLES =";
  for (unsigned d=0; d<dim; ++d)
    f << " \"" << m.vn[d] << "\"";
  f << endl;


  // check volume mesh
  const unsigned Nelem = m.e();
  for (unsigned c=0; c<Nelem; ++c) {
    if (!(c%100000))
      cout << "geometry checking progress: " << (int)(100.*((double)c)/((double)Nelem)) << "%" << endl;
    if (m.vz[0].e2n[c].n.size()!=dim+1)
      continue;

    // get positive/negative volumes
    cellgeom(m,c,&No_loc[0],&vol,&inc_min);
    n_positivesize += (vol>0.? 1:0);
    n_negativesize += (vol<0.? 1:0);

    // output elements with normals (2 zones per element)
    if (vol<=0.) {
      std::vector< unsigned >& n = m.vz[0].e2n[c].n;
      if (dim>2) {
        f << "ZONE T=\"Element_" << c << "\", N=4, E=1, ZONETYPE=FETETRAHEDRON, DATAPACKING=POINT" << endl;
        for (unsigned i=0; i<4; ++i)
          f << m.vv[0][n[i]] << " " << m.vv[1][n[i]] << " " << m.vv[2][n[i]] << endl;
        f << "1 2 3 4" << endl;
      }
      else {
        f << "ZONE T=\"Element_" << c << "\", N=3, E=1, ZONETYPE=FETRIANGLE, DATAPACKING=POINT" << endl;
        for (unsigned i=0; i<3; ++i)
          f << m.vv[0][n[i]] << " " << m.vv[1][n[i]] << endl;
        f << "1 2 3" << endl;
      }
    }

  }
  if (Nelem%100000)
    cout << "geometry checking progress: 100%" << endl;
  cout << "geometry element volumes summary (+/-): " << n_positivesize << "/" << n_negativesize << " of " << Nelem << endl;


  // check boundaries mesh (output all face normals)

  // create normals beggining/ending points and normalize
  // lengths because tecplot can't see small lengths
  const double nlength = 1.e-1;  // intended normals length
  vector< vector< double > > bxi(m.vz.size()), byi(m.vz.size()), bzi(m.vz.size()),
                             bxf(m.vz.size()), byf(m.vz.size()), bzf(m.vz.size());
  const vector< double >& vx = m.vv[0];
  const vector< double >& vy = m.vv[1];
  for (unsigned t=1; t<m.vz.size(); ++t) {
    vector< melem >& vbe = m.vz[t].e2n;     // reference
    for (unsigned c=0; c<vbe.size(); ++c) {
      if (dim==2) {
        if (vbe[c].n.size()==2) {
          const unsigned n1 = vbe[c].n[0];
          const unsigned n2 = vbe[c].n[1];
          bxi[t].push_back((vx[n1]+vx[n2])/2.);
          byi[t].push_back((vy[n1]+vy[n2])/2.);
          bxf[t].push_back(vy[n1]-vy[n2]);  // normal = (0,0,1) x (vx,vy,0)
          byf[t].push_back(vx[n2]-vx[n1]);
        }
      }
      else if (dim==3) {
        const vector< double >& vz = m.vv[2];
        if (vbe[c].n.size()==3) {
          const unsigned n1 = vbe[c].n[0];
          const unsigned n2 = vbe[c].n[1];
          const unsigned n3 = vbe[c].n[2];
          bxi[t].push_back((vx[n1]+vx[n2]+vx[n3])/3.);
          byi[t].push_back((vy[n1]+vy[n2]+vy[n3])/3.);
          bzi[t].push_back((vz[n1]+vz[n2]+vz[n3])/3.);
          bxf[t].push_back((vy[n2]-vy[n1])*(vz[n3]-vz[n1])-(vz[n2]-vz[n1])*(vy[n3]-vy[n1]));
          byf[t].push_back((vz[n2]-vz[n1])*(vx[n3]-vx[n1])-(vx[n2]-vx[n1])*(vz[n3]-vz[n1]));
          bzf[t].push_back((vx[n2]-vx[n1])*(vy[n3]-vy[n1])-(vy[n2]-vy[n1])*(vx[n3]-vx[n1]));
        }
        else if (vbe[c].n.size()==4) {
          const unsigned n1 = vbe[c].n[0];
          const unsigned n2 = vbe[c].n[1];
          const unsigned n3 = vbe[c].n[2];
          const unsigned n4 = vbe[c].n[3];
          bxi[t].push_back((vx[n1]+vx[n2]+vx[n3])/3.);
          byi[t].push_back((vy[n1]+vy[n2]+vy[n3])/3.);
          bzi[t].push_back((vz[n1]+vz[n2]+vz[n3])/3.);
          bxf[t].push_back((vy[n2]-vy[n1])*(vz[n3]-vz[n1])-(vz[n2]-vz[n1])*(vy[n3]-vy[n1]));
          byf[t].push_back((vz[n2]-vz[n1])*(vx[n3]-vx[n1])-(vx[n2]-vx[n1])*(vz[n3]-vz[n1]));
          bzf[t].push_back((vx[n2]-vx[n1])*(vy[n3]-vy[n1])-(vy[n2]-vy[n1])*(vx[n3]-vx[n1]));
          bxi[t].push_back((vx[n1]+vx[n3]+vx[n4])/3.);
          byi[t].push_back((vy[n1]+vy[n3]+vy[n4])/3.);
          bzi[t].push_back((vz[n1]+vz[n3]+vz[n4])/3.);
          bxf[t].push_back((vy[n3]-vy[n1])*(vz[n4]-vz[n1])-(vz[n3]-vz[n1])*(vy[n4]-vy[n1]));
          byf[t].push_back((vz[n3]-vz[n1])*(vx[n4]-vx[n1])-(vx[n3]-vx[n1])*(vz[n4]-vz[n1]));
          bzf[t].push_back((vx[n3]-vx[n1])*(vy[n4]-vy[n1])-(vy[n3]-vy[n1])*(vx[n4]-vx[n1]));
        }
      }
    }
  }

  // normalize normals' length (normal is in b[xyz]f) and set final position
  for (unsigned t=1; t<m.vz.size(); ++t) {
    for (unsigned i=0; i<bxi[t].size(); ++i)
      if (dim==2) {
        const double length = sqrt(bxf[t][i]*bxf[t][i]+byf[t][i]*byf[t][i]);
        bxf[t][i] = bxi[t][i] + bxf[t][i]*nlength/length;
        byf[t][i] = byi[t][i] + byf[t][i]*nlength/length;
      }
      else if (dim==3) {
        const double length = sqrt(bxf[t][i]*bxf[t][i]+byf[t][i]*byf[t][i]+bzf[t][i]*bzf[t][i]);
        bxf[t][i] = bxi[t][i] + bxf[t][i]*nlength/length;
        byf[t][i] = byi[t][i] + byf[t][i]*nlength/length;
        bzf[t][i] = bzi[t][i] + bzf[t][i]*nlength/length;
      }
  }
  for (unsigned t=1; t<m.vz.size(); ++t) {
    f << "ZONE T=\"" << m.vz[t].n << "\", N=" << 2*bxi[t].size() << ", E=" << bxi[t].size() << ", ZONETYPE=FELINESEG, DATAPACKING=POINT" << endl;
    if (dim==2)
      for (unsigned i=0; i<bxi[t].size(); ++i)
        f << bxi[t][i] << " " << byi[t][i] << endl
          << bxf[t][i] << " " << byf[t][i] << endl;
    else if (dim==3)
      for (unsigned i=0; i<bxi[t].size(); ++i)
        f << bxi[t][i] << " " << byi[t][i] << " " << bzi[t][i] << endl
          << bxf[t][i] << " " << byf[t][i] << " " << bzf[t][i] << endl;
    for (unsigned i=0; i<2*bxi[t].size(); i+=2)
      f << 1+i << " " << 2+i << endl;
  }


  // finished
  f.close();
  cout << "geometry checking finished" << endl;
}


// if point is inside triangle
bool t_geo::ispointintriangle(
  const double x,  const double y,
  const double x1, const double y1,
  const double x2, const double y2,
  const double x3, const double y3 )
{
  // vectors
  double v0[2] = { x3-x1, y3-y1 };
  double v1[2] = { x2-x1, y2-y1 };
  double v2[2] = { x -x1, y -y1 };

  // dot products
  const double dot00 = v0[0]*v0[0] + v0[1]*v0[1];
  const double dot01 = v0[0]*v1[0] + v0[1]*v1[1];
  const double dot02 = v0[0]*v2[0] + v0[1]*v2[1];
  const double dot11 = v1[0]*v1[0] + v1[1]*v1[1];
  const double dot12 = v1[0]*v2[0] + v1[1]*v2[1];

  // barycentric coordinates
  const double invdenom = 1. / (dot00*dot11 - dot01*dot01);
  const double u = (dot11*dot02 - dot01*dot12)*invdenom;
  const double v = (dot00*dot12 - dot01*dot02)*invdenom;

  // check if point is in triangle
  return ((u>=0) && (v>=0) && (u+v<1+eps));
}


// if point is inside tetrahedron
bool t_geo::ispointintetrahedron(
  const double x,  const double y,  const double z,
  const double x1, const double y1, const double z1,
  const double x2, const double y2, const double z2,
  const double x3, const double y3, const double z3,
  const double x4, const double y4, const double z4 )
{
  // notes:
  //  - barycentric coordinates are bi = Di/D0, i:[1,4]
  //  - condition must hold: sum(bi)=1, or sum(Di)=D0
  //  - if bi=0, then P lies on boundary i (the one formed by the points but Vi)
  //  - if bi>0, then P is inside boundary i (same side as Vi), and vice-versa
  //  - if P is inside all 4 boundaries, then it is inside tetrahedron
  //  - this can be extended to simplexes of any dimension

  //      |x1 y1 z1 1|       |x  y  z  1|       |x1 y1 z1 1|       |x1 y1 z1 1|       |x1 y1 z1 1|
  // D0 = |x2 y2 z2 1|  D1 = |x2 y2 z2 1|  D2 = |x  y  z  1|  D3 = |x2 y2 z2 1|  D4 = |x2 y2 z2 1|
  //      |x3 y3 z3 1|       |x3 y3 z3 1|       |x3 y3 z3 1|       |x  y  z  1|       |x3 y3 z3 1|
  //      |x4 y4 z4 1|       |x4 y4 z4 1|       |x4 y4 z4 1|       |x4 y4 z4 1|       |x  y  z  1|
  const double D0 =
    - (x2*(y3*z4-z3*y4)-y2*(x3*z4-z3*x4)+z2*(x3*y4-y3*x4))
    + (x1*(y3*z4-z3*y4)-y1*(x3*z4-z3*x4)+z1*(x3*y4-y3*x4))
    - (x1*(y2*z4-z2*y4)-y1*(x2*z4-z2*x4)+z1*(x2*y4-y2*x4))
    + (x1*(y2*z3-z2*y3)-y1*(x2*z3-z2*x3)+z1*(x2*y3-y2*x3));

  // check if tetrahedron is degenerate (coplanar points, zero volume)
  if (abs(D0)<eps) {
    cout << "tetrahedron is degenerate: D0=" << D0 << endl;
    return false;
  }

  const double D1 =
    - (x2*(y3*z4-z3*y4)-y2*(x3*z4-z3*x4)+z2*(x3*y4-y3*x4))
    + (x *(y3*z4-z3*y4)-y *(x3*z4-z3*x4)+z *(x3*y4-y3*x4))
    - (x *(y2*z4-z2*y4)-y *(x2*z4-z2*x4)+z *(x2*y4-y2*x4))
    + (x *(y2*z3-z2*y3)-y *(x2*z3-z2*x3)+z *(x2*y3-y2*x3));
  const double D2 =
    - (x *(y3*z4-z3*y4)-y *(x3*z4-z3*x4)+z *(x3*y4-y3*x4))
    + (x1*(y3*z4-z3*y4)-y1*(x3*z4-z3*x4)+z1*(x3*y4-y3*x4))
    - (x1*(y *z4-z *y4)-y1*(x *z4-z *x4)+z1*(x *y4-y *x4))
    + (x1*(y *z3-z *y3)-y1*(x *z3-z *x3)+z1*(x *y3-y *x3));
  const double D3 =
    - (x2*(y *z4-z *y4)-y2*(x *z4-z *x4)+z2*(x *y4-y *x4))
    + (x1*(y *z4-z *y4)-y1*(x *z4-z *x4)+z1*(x *y4-y *x4))
    - (x1*(y2*z4-z2*y4)-y1*(x2*z4-z2*x4)+z1*(x2*y4-y2*x4))
    + (x1*(y2*z -z2*y )-y1*(x2*z -z2*x )+z1*(x2*y -y2*x ));
  const double D4 =
    - (x2*(y3*z -z3*y )-y2*(x3*z -z3*x )+z2*(x3*y -y3*x ))
    + (x1*(y3*z -z3*y )-y1*(x3*z -z3*x )+z1*(x3*y -y3*x ))
    - (x1*(y2*z -z2*y )-y1*(x2*z -z2*x )+z1*(x2*y -y2*x ))
    + (x1*(y2*z3-z2*y3)-y1*(x2*z3-z2*x3)+z1*(x2*y3-y2*x3));

  // check if barycentric coordinates are correct
  if (abs(D1+D2+D3+D4 - D0)>eps) {
    cout << "tetrahedron barycentric coordinates not correct:" <<
    " b1=" << D1/D0 << " b2=" << D2/D0 << " b3=" << D3/D0 << " b4=" << D4/D0 << endl;
    return false;
  }

  // if any barycentric coordinate is negative, point is outside
  if (D0>0.) {
    if (D1>0. && D2>0. && D3>0. && D4>0.)
      return true;
  }
  else {
    if (D1<0. && D2<0. && D3<0. && D4<0.)
      return true;
  }
  return false;
}


bool t_geo::cellgeom_ok(const vector< vector< double > >& coords)
{
  const int Nvtcell = (int) coords.size();
  const int Ndim    = (int) (Nvtcell? coords[0].size():0);

  #define cf_assert_desc( m , c ) if (!(c)) return false;

  double d03[3],d13[3],d23[3],d01[3],d21[3],d31[3];
  vector< vector< double > > norm(Nvtcell,vector< double >(Ndim,0.));

  if (Ndim==2 && Nvtcell==3) {
    // triangles

    for (int j0=0; j0<Nvtcell; ++j0) {
      const int j1 = (j0+1)%Nvtcell;
      const int j2 = (j0+2)%Nvtcell;

      // inwards-facing normals (and vectors j2->j0)
      norm[j0][0] = coords[j1][1] - coords[j2][1];
      norm[j0][1] = coords[j2][0] - coords[j1][0];
      const double rvec0 = coords[j0][0] - coords[j2][0];
      const double rvec1 = coords[j0][1] - coords[j2][1];
      cf_assert_desc( "normal not inward-facing!",
        norm[j0][0]*rvec0 + norm[j0][1]*rvec1 >= 0. );
    }

  }
  else if (Ndim==3 && Nvtcell==4) {
    // tetrahedra

    // inward-facing normals
    for (int id=0; id<Ndim; ++id) {
      d23[id] = coords[2][id] - coords[3][id];
      d03[id] = coords[0][id] - coords[3][id];
      d13[id] = coords[1][id] - coords[3][id];
      d01[id] = coords[0][id] - coords[1][id];
      d21[id] = coords[2][id] - coords[1][id];
      d31[id] = -d13[id];
    }
#define vecprd3(v1,v2,v) { v[0] = v1[1]*v2[2]-v1[2]*v2[1]; v[1] = v1[2]*v2[0]-v1[0]*v2[2]; v[2] = v1[0]*v2[1]-v1[1]*v2[0]; }
    vecprd3(d23,d13,norm[0]);
    vecprd3(d03,d23,norm[1]);
    vecprd3(d01,d31,norm[2]);
    vecprd3(d21,d01,norm[3]);
#undef vecprd3
    cf_assert_desc( "normal not inward-facing!",
      d03[0] * norm[0][0] +
      d03[1] * norm[0][1] +
      d03[2] * norm[0][2] >= 0.);

    // halve magnitude
    for (int inc=0; inc<Nvtcell; ++inc)
      for (int id=0; id<Ndim; ++id)
        norm[inc][id] *= .5;

  }
  else {
    return true;
  }

  // volume
  double vol = 0.;
  for (int inc=0; inc<Nvtcell; ++inc)
    for (int id=0; id<Ndim; ++id)
      vol += coords[inc][id] * norm[inc][id];
  cf_assert_desc("negative volume on element!",vol>0.);
  #undef cf_assert_desc

  return true;
}


void t_geo::cellgeom(
  mmesh& m,
  int ic,
  struct local_node_struct *No_loc,
  double *vol,
  int *inc_min )
{
  const int Ndim    = (int) m.d();
  const int Nvtcell = Ndim+1;
  const double dNvtfce = (double) Ndim;
  #define cf_assert_desc( m , c ) if (!(c)) { std::cerr << "element " << ic << ": " << (m) << std::endl; }

  std::vector< std::vector< double > > coords(Nvtcell,std::vector< double >(Ndim,0.));

  double d03[3],d13[3],d23[3],d01[3],d21[3],d31[3];

  if (m.vz[0].e2n[ic].n.size() > (unsigned) Nvtcell)
    return;
  for (int inc=0; inc<Nvtcell; ++inc) {
    const int n = m.vz[0].e2n[ic].n[inc];
    No_loc[inc].node = n;
    coords[inc][0] = m.vv[0][n];
    coords[inc][1] = m.vv[1][n];
    if (Ndim==3)
      coords[inc][2] = m.vv[2][n];
    for (int id=0; id<Ndim; ++id)
      No_loc[inc].norm[id] = 0.;
  }
  *inc_min = 0;

  if (Ndim==2) {
    // triangles

    for (int j0=0; j0<Nvtcell; ++j0) {
      const int j1 = (j0+1)%Nvtcell;
      const int j2 = (j0+2)%Nvtcell;

      // inwards-facing normals (and vectors j2->j0)
      No_loc[j0].norm[0] = coords[j1][1] - coords[j2][1];
      No_loc[j0].norm[1] = coords[j2][0] - coords[j1][0];
      const double rvec0 = coords[j0][0] - coords[j2][0];
      const double rvec1 = coords[j0][1] - coords[j2][1];
      cf_assert_desc( "normal not inward-facing!",
        No_loc[j0].norm[0]*rvec0 + No_loc[j0].norm[1]*rvec1 >= 0. );

      // set square of magnitude and face with smallest area
      No_loc[j0].norm2 =
        No_loc[j0].norm[0] * No_loc[j0].norm[0] +
        No_loc[j0].norm[1] * No_loc[j0].norm[1];
      if (No_loc[j0].norm2 < No_loc[*inc_min].norm2)
        *inc_min=j0;
    }

  }
  else {
    // tetrahedra

    // inward-facing normals
    for (int id=0; id<Ndim; ++id) {
      d23[id] = coords[2][id] - coords[3][id];
      d03[id] = coords[0][id] - coords[3][id];
      d13[id] = coords[1][id] - coords[3][id];
      d01[id] = coords[0][id] - coords[1][id];
      d21[id] = coords[2][id] - coords[1][id];
      d31[id] = -d13[id];
    }
    vecprd3(d23,d13,No_loc[0].norm);
    vecprd3(d03,d23,No_loc[1].norm);
    vecprd3(d01,d31,No_loc[2].norm);
    vecprd3(d21,d01,No_loc[3].norm);
    cf_assert_desc( "normal not inward-facing!",
      d03[0] * No_loc[0].norm[0] +
      d03[1] * No_loc[0].norm[1] +
      d03[2] * No_loc[0].norm[2] >= 0.);

    // halve magnitude
    for (int inc=0; inc<Nvtcell; ++inc)
      for (int id=0; id<Ndim; ++id)
        No_loc[inc].norm[id] *= .5;

    // set square of magnitude and face with smallest area
    for (int inc=0; inc<Nvtcell; ++inc) {
      No_loc[inc].norm2 =
        No_loc[inc].norm[0]*No_loc[inc].norm[0] +
        No_loc[inc].norm[1]*No_loc[inc].norm[1] +
        No_loc[inc].norm[2]*No_loc[inc].norm[2];
      if (No_loc[inc].norm2 < No_loc[*inc_min].norm2)
        *inc_min = inc;
    }

  }

  // volume
  *vol = 0.;
  for (int inc=0; inc<Nvtcell; ++inc)
    for (int id=0; id<Ndim; ++id)
      *vol += coords[inc][id] * No_loc[inc].norm[id];
  *vol = *vol/((double)Ndim * dNvtfce);
  cf_assert_desc("negative volume on element!",*vol>0.);
  #undef cf_assert_desc
}


void t_geo::vecprd3(double *v1, double *v2, double *v)
{
  v[0] = v1[1]*v2[2]-v1[2]*v2[1];
  v[1] = v1[2]*v2[0]-v1[0]*v2[2];
  v[2] = v1[0]*v2[1]-v1[1]*v2[0];
}

