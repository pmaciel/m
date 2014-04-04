
#include "pillaz.h"



double pillaz::pillaz_CalcSizeTriangle(const unsigned n1, const unsigned n2, const unsigned n3)
{
  double c[3][2];
  getQuantityVector(COORD,n1,2,c[0]);
  getQuantityVector(COORD,n2,2,c[1]);
  getQuantityVector(COORD,n3,2,c[2]);
  return c[0][0]*(c[1][1]-c[2][1])
       + c[1][0]*(c[2][1]-c[0][1])
       + c[2][0]*(c[0][1]-c[1][1]);
}


double pillaz::pillaz_CalcSizeTetrahedron(const unsigned n1, const unsigned n2, const unsigned n3, const unsigned n4)
{
  double c[4][3];
  getQuantityVector(COORD,n1,3,c[0]);
  getQuantityVector(COORD,n2,3,c[1]);
  getQuantityVector(COORD,n3,3,c[2]);
  getQuantityVector(COORD,n4,3,c[3]);

  double v1[3], v2[3], v3[3], v4[3];
  for (int d=0; d<3; ++d) {
    v1[d] = c[1][d] - c[0][d];
    v2[d] = c[2][d] - c[0][d];
    v3[d] = c[3][d] - c[0][d];
  }
  pillaz_CalcCrossProduct_3D(v4,v1,v2);
  return (pillaz_CalcVectorScalarProduct(3,v3,v4)/6.);
}


double pillaz::pillaz_CalcElmSize(const unsigned iz, const unsigned ie)
{
  pillaz_elmtype_t et(getElmType(iz,ie));
  std::vector< int > n(pillaz_getElmNNodes(et),-1);
  getElmNodes(iz,ie,&n[0]);

  switch (et) {
  case ELM_TRIANGLE:    return pillaz_CalcSizeTriangle(n[0],n[1],n[2]);
  case ELM_TETRAHEDRON: return pillaz_CalcSizeTetrahedron(n[0],n[1],n[2],n[3]);
  case ELM_WEDGE:       return pillaz_CalcSizeTetrahedron(n[0],n[1],n[2],n[5])
                             + pillaz_CalcSizeTetrahedron(n[0],n[1],n[5],n[4])
                             + pillaz_CalcSizeTetrahedron(n[0],n[4],n[5],n[3]);
  case ELM_QUAD:        return pillaz_CalcSizeTriangle(n[0],n[1],n[2])
                             + pillaz_CalcSizeTriangle(n[0],n[2],n[3]);
  case ELM_BRICK:       return pillaz_CalcSizeTetrahedron(n[0],n[1],n[3],n[5])
                             + pillaz_CalcSizeTetrahedron(n[0],n[3],n[6],n[5])
                             + pillaz_CalcSizeTetrahedron(n[0],n[3],n[2],n[6])
                             + pillaz_CalcSizeTetrahedron(n[0],n[5],n[6],n[4])
                             + pillaz_CalcSizeTetrahedron(n[3],n[6],n[5],n[7]);
  case ELM_PYRAMID:     return pillaz_CalcSizeTetrahedron(n[0],n[1],n[3],n[4])
                             + pillaz_CalcSizeTetrahedron(n[0],n[3],n[2],n[4]);
  case ELM_EDGE:
  case ELM_UNDEFINED:
  case ALL_ELEMENTS:
  default: break;
  }
  return 0.;
}


void pillaz::pillaz_RandomElmPosition(const int zone, const int elm, double *p)
{
  // random numbers
  const double
    rand1 = pillaz_RandomDouble(),
    rand2 = pillaz_RandomDouble(),
    rand3 = pillaz_RandomDouble(),
    rand4 = pillaz_RandomDouble();

  // get element nodes and coordinates
  std::vector< int > n(pillaz_getElmNNodes(getElmType(zone,elm)),-1);
  std::vector< std::vector< double > > c(n.size(),std::vector< double >(fp.numDim,0.));
  getElmNodes(zone,elm,&n[0]);
  for (size_t i=0; i<n.size(); ++i)
    getQuantityVector(COORD,n[i],fp.numDim,&(c[i])[0]);

  switch (getElmType(zone,elm)) {
  case ELM_TRIANGLE: {

    for (int d=0; d<fp.numDim; ++d)
      p[d] = c[0][d] + rand1*(c[1][d] - c[0][d]) + (1.-rand1)*rand2*(c[2][d] - c[0][d]);

  } break;
  case ELM_QUAD: {

    for (int d=0; d<fp.numDim; ++d)
      p[d] = rand3 * ( c[0][d] + rand1*(c[1][d] - c[0][d])) + (1.-rand3) * ( c[2][d] + rand2*(c[3][d] - c[2][d]));

  } break;
  case ELM_TETRAHEDRON: {

    for (int d=0; d<fp.numDim; ++d)
      p[d] = c[0][d] + rand1*(c[1][d] - c[0][d]) + (1.-rand1)*rand2*(c[2][d] - c[0][d]) + (1.-rand1-(1.-rand1)*rand2)*rand3*(c[3][d] - c[0][d]);

  } break;
  case ELM_BRICK: {

    for (int d=0; d<fp.numDim; ++d)
      p[d] = rand4*(rand3*(c[0][d] + rand1*(c[4][d] - c[0][d])) + (1.-rand3)*(c[1][d] + rand2*(c[5][d] - c[1][d]))) +(1.-rand4)*(rand3*(c[2][d] + rand1*(c[6][d] - c[2][d])) + (1.-rand3)*(c[3][d] + rand2*(c[7][d] - c[3][d])));

  } break;
  case ELM_WEDGE: {

    for (int d=0; d<fp.numDim; ++d)
      p[d] = rand3*c[0][d] + rand3*rand1*(c[1][d] - c[0][d]) + rand3*(1.-rand1)*rand2*(c[2][d] - c[0][d]) + (1.-rand3)*c[3][d] + (1.-rand3)*rand1*(c[4][d] - c[3][d]) + (1.-rand3)*(1.-rand1)*rand2*(c[5][d] - c[3][d]);

  } break;
  case ELM_PYRAMID: {

    for (int d=0; d<fp.numDim; ++d)
      p[d] = rand4*c[4][d] + (1.-rand4)*(rand3*(c[0][d] + rand1*(c[1][d] - c[0][d])) + (1.-rand3)*(c[2][d] + rand2*(c[3][d] - c[2][d])));

  } break;
  default: break;
  }

}


void pillaz::pillaz_getElmFaceMiddlePoint(const int iz, const int ie, const int iface, double *fmp)
{
  // initialize
  pillaz_elmtype_t ft = pillaz_getElmFaceType(getElmType(iz,ie),iface);
  std::vector< int > fn(pillaz_getElmNNodes(ft),-1);
  std::vector< std::vector< double > > c(fn.size(),std::vector< double >(fp.numDim,0.));

  // get element nodes and coordinates
  pillaz_getElmFaceNodes(iz,ie,iface,&fn[0]);
  for (size_t i=0; i<fn.size(); ++i)
    getQuantityVector(COORD,fn[i],fp.numDim,&(c[i])[0]);

  // calculate middle point
  for (int d=0; d<fp.numDim; ++d) {
    fmp[d] = 0.;
    for (size_t i=0; i<fn.size(); ++i)
      fmp[d] += c[i][d];
    fmp[d] /= (double) fn.size();
  }
}


pillaz_elmtype_t pillaz::pillaz_getElmFaceType(const pillaz_elmtype_t et, const int iface)
{
  return (et==ELM_TRIANGLE?    ELM_EDGE     :
         (et==ELM_TETRAHEDRON? ELM_TRIANGLE :
         (et==ELM_QUAD?        ELM_EDGE     :
         (et==ELM_BRICK?       ELM_QUAD     :
         (et==ELM_WEDGE?       (iface<3? ELM_QUAD : ELM_TRIANGLE) :
         (et==ELM_PYRAMID?     (iface<1? ELM_QUAD : ELM_TRIANGLE) :
                               ELM_UNDEFINED ))))));
}


void pillaz::pillaz_getElmFaceNodes(const int iz, const int ie, const int iface, int *fnodes)
{
  // get element nodes
  const pillaz_elmtype_t et = getElmType(iz,ie);
  std::vector< int > en(pillaz_getElmNNodes(et),-1);
  getElmNodes(iz,ie,&en[0]);

  // get face nodes
  switch (et) {
  case ELM_TRIANGLE:
    switch (iface) {
    case 0: { fnodes[0]=en[0]; fnodes[1]=en[1]; } break;
    case 1: { fnodes[0]=en[1]; fnodes[1]=en[2]; } break;
    case 2: { fnodes[0]=en[2]; fnodes[1]=en[0]; } break;
    } break;
  case ELM_TETRAHEDRON:
    switch (iface) {
    case 0: { fnodes[0]=en[1]; fnodes[1]=en[0]; fnodes[2]=en[2]; } break;
    case 1: { fnodes[0]=en[0]; fnodes[1]=en[1]; fnodes[2]=en[3]; } break;
    case 2: { fnodes[0]=en[1]; fnodes[1]=en[2]; fnodes[2]=en[3]; } break;
    case 3: { fnodes[0]=en[2]; fnodes[1]=en[0]; fnodes[2]=en[3]; } break;
    } break;
  case ELM_QUAD:
    switch (iface) {
    case 0: { fnodes[0]=en[0]; fnodes[1]=en[1]; } break;
    case 1: { fnodes[0]=en[1]; fnodes[1]=en[2]; } break;
    case 2: { fnodes[0]=en[2]; fnodes[1]=en[3]; } break;
    case 3: { fnodes[0]=en[3]; fnodes[1]=en[0]; } break;
    } break;
  case ELM_BRICK:
    switch (iface) {
    case 0: { fnodes[0]=en[0]; fnodes[1]=en[1]; fnodes[2]=en[5]; fnodes[3]=en[4]; } break;
    case 1: { fnodes[0]=en[1]; fnodes[1]=en[3]; fnodes[2]=en[7]; fnodes[3]=en[5]; } break;
    case 2: { fnodes[0]=en[3]; fnodes[1]=en[2]; fnodes[2]=en[6]; fnodes[3]=en[7]; } break;
    case 3: { fnodes[0]=en[2]; fnodes[1]=en[0]; fnodes[2]=en[4]; fnodes[3]=en[6]; } break;
    case 4: { fnodes[0]=en[1]; fnodes[1]=en[0]; fnodes[2]=en[2]; fnodes[3]=en[3]; } break;
    case 5: { fnodes[0]=en[4]; fnodes[1]=en[5]; fnodes[2]=en[7]; fnodes[3]=en[6]; } break;
    } break;
  case ELM_WEDGE:
    switch (iface) {
    case 0: { fnodes[0]=en[0]; fnodes[1]=en[1]; fnodes[2]=en[4]; fnodes[3]=en[3]; } break;
    case 1: { fnodes[0]=en[1]; fnodes[1]=en[2]; fnodes[2]=en[5]; fnodes[3]=en[4]; } break;
    case 2: { fnodes[0]=en[2]; fnodes[1]=en[0]; fnodes[2]=en[3]; fnodes[3]=en[5]; } break;
    case 3: { fnodes[0]=en[0]; fnodes[1]=en[2]; fnodes[2]=en[1]; } break;
    case 4: { fnodes[0]=en[3]; fnodes[1]=en[4]; fnodes[2]=en[5]; } break;
    } break;
  case ELM_PYRAMID:
    switch (iface) {
    case 0: { fnodes[0]=en[0]; fnodes[1]=en[2]; fnodes[2]=en[3]; fnodes[3]=en[1]; } break;
    case 1: { fnodes[0]=en[0]; fnodes[1]=en[1]; fnodes[2]=en[4]; } break;
    case 2: { fnodes[0]=en[1]; fnodes[1]=en[3]; fnodes[2]=en[4]; } break;
    case 3: { fnodes[0]=en[3]; fnodes[1]=en[2]; fnodes[2]=en[4]; } break;
    case 4: { fnodes[0]=en[2]; fnodes[1]=en[0]; fnodes[2]=en[4]; } break;
    } break;
  case ALL_ELEMENTS:
  default: break;
  }
}


int pillaz::pillaz_getElmNFaces(const pillaz_elmtype_t et)
{
  return (et==ELM_TRIANGLE?    3 :
         (et==ELM_TETRAHEDRON? 4 :
         (et==ELM_WEDGE?       5 :
         (et==ELM_QUAD?        4 :
         (et==ELM_BRICK?       6 :
         (et==ELM_PYRAMID?     5 : 0 ))))));
}


int pillaz::pillaz_getElmNNodes(const pillaz_elmtype_t et)
{
  return (et==ELM_TRIANGLE?    3 :
         (et==ELM_TETRAHEDRON? 4 :
         (et==ELM_WEDGE?       6 :
         (et==ELM_QUAD?        4 :
         (et==ELM_BRICK?       8 :
         (et==ELM_PYRAMID?     5 : 0 ))))));
}


void pillaz::pillaz_CalcElmFaceNormal(const int zone, const int elm, const int eface, double *normal)
{
  // initialize
  for (int d=0; d<fp.numDim; ++d)
    normal[d] = 0.;
  const pillaz_elmtype_t
    et = getElmType(zone,elm),
    ft = pillaz_getElmFaceType(et,eface);
  std::vector< int > fn(pillaz_getElmNNodes(ft),-1);
  std::vector< std::vector< double > > c(fn.size(),std::vector< double >(fp.numDim,0.));

  // get element nodes and coordinates
  pillaz_getElmFaceNodes(zone,elm,eface,&fn[0]);
  for (size_t i=0; i<fn.size(); ++i)
    getQuantityVector(COORD,fn[i],fp.numDim,&(c[i])[0]);

  // calculate element face normal
  if (et==ELM_TRIANGLE) {
    normal[0] = c[0][1] - c[1][1];
    normal[1] = c[1][0] - c[0][0];
  }
  else if ((et==ELM_TETRAHEDRON)
        || (et==ELM_WEDGE   && eface>2)
        || (et==ELM_PYRAMID && eface>0)) {
    normal[0] = 0.5*((c[2][1] - c[0][1]) * (c[1][2] - c[0][2]) - (c[1][1] - c[0][1]) * (c[2][2] - c[0][2]));
    normal[1] = 0.5*((c[2][2] - c[0][2]) * (c[1][0] - c[0][0]) - (c[1][2] - c[0][2]) * (c[2][0] - c[0][0]));
    normal[2] = 0.5*((c[2][0] - c[0][0]) * (c[1][1] - c[0][1]) - (c[1][0] - c[0][0]) * (c[2][1] - c[0][1]));
  }
  else if(et==ELM_QUAD) {
    normal[0] = 2.*(c[0][1] - c[1][1]);
    normal[1] = 2.*(c[1][0] - c[0][0]);
  }
  else if (et==ELM_BRICK
       || (et==ELM_WEDGE   && eface<=2)
       || (et==ELM_PYRAMID && eface==0)) {
    normal[0] = ((c[3][1] - c[0][1]) * (c[1][2] - c[0][2]) -(c[1][1] - c[0][1]) * (c[3][2] - c[0][2]));
    normal[1] = ((c[3][2] - c[0][2]) * (c[1][0] - c[0][0]) -(c[1][2] - c[0][2]) * (c[3][0] - c[0][0]));
    normal[2] = ((c[3][0] - c[0][0]) * (c[1][1] - c[0][1]) -(c[1][0] - c[0][0]) * (c[3][1] - c[0][1]));
  }

}

