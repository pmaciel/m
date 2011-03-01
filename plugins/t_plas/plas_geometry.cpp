
#include "plas.h"



double plas::plas_CalcSizeTriangle(const unsigned n1, const unsigned n2, const unsigned n3)
{
  double c[3][2];
  c[0][0]=plas_getQuantity(COORD_X,n1);  c[0][1]=plas_getQuantity(COORD_Y,n1);
  c[1][0]=plas_getQuantity(COORD_X,n2);  c[1][1]=plas_getQuantity(COORD_Y,n2);
  c[2][0]=plas_getQuantity(COORD_X,n3);  c[2][1]=plas_getQuantity(COORD_Y,n3);
  return c[0][0]*(c[1][1]-c[2][1])
       + c[1][0]*(c[2][1]-c[0][1])
       + c[2][0]*(c[0][1]-c[1][1]);
}


double plas::plas_CalcSizeTetrahedron(const unsigned n1, const unsigned n2, const unsigned n3, const unsigned n4)
{
  double p1[3], v1[3], v2[3], v3[3], v4[3];
  p1[0] = plas_getQuantity(COORD_X,n1);
  p1[1] = plas_getQuantity(COORD_Y,n1);
  p1[2] = plas_getQuantity(COORD_Z,n1);
  v1[0] = plas_getQuantity(COORD_X,n2) - p1[0];
  v1[1] = plas_getQuantity(COORD_Y,n2) - p1[1];
  v1[2] = plas_getQuantity(COORD_Z,n2) - p1[2];
  v2[0] = plas_getQuantity(COORD_X,n3) - p1[0];
  v2[1] = plas_getQuantity(COORD_Y,n3) - p1[1];
  v2[2] = plas_getQuantity(COORD_Z,n3) - p1[2];
  v3[0] = plas_getQuantity(COORD_X,n4) - p1[0];
  v3[1] = plas_getQuantity(COORD_Y,n4) - p1[1];
  v3[2] = plas_getQuantity(COORD_Z,n4) - p1[2];
  plas_CalcCrossProduct_3D(v4,v1,v2);
  return (plas_CalcVectorScalarProduct(3,v3,v4)/6.);
}


double plas::plas_CalcSizePrism(const unsigned n1, const unsigned n2, const unsigned n3, const unsigned n4, const unsigned n5, const unsigned n6)
{
  return plas_CalcSizeTetrahedron(n1,n2,n3,n6)
       + plas_CalcSizeTetrahedron(n1,n2,n6,n5)
       + plas_CalcSizeTetrahedron(n1,n5,n6,n4);
}


double plas::plas_CalcSizeQuadrilateral(const unsigned n1, const unsigned n2, const unsigned n3, const unsigned n4)
{
  return plas_CalcSizeTriangle(n1,n2,n3)
       + plas_CalcSizeTriangle(n1,n3,n4);
}


double plas::plas_CalcSizeHexahedron(const unsigned n1, const unsigned n2, const unsigned n3, const unsigned n4, const unsigned n5, const unsigned n6, const unsigned n7, const unsigned n8)
{
  return plas_CalcSizeTetrahedron(n1,n2,n4,n6)
       + plas_CalcSizeTetrahedron(n1,n4,n7,n6)
       + plas_CalcSizeTetrahedron(n1,n4,n3,n7)
       + plas_CalcSizeTetrahedron(n1,n6,n7,n5)
       + plas_CalcSizeTetrahedron(n4,n7,n6,n8);
}


double plas::plas_CalcSizePyramid(const unsigned n1, const unsigned n2, const unsigned n3, const unsigned n4, const unsigned n5)
{
  return plas_CalcSizeTetrahedron(n1,n2,n4,n5)
       + plas_CalcSizeTetrahedron(n1,n4,n3,n5);
}

void plas::plas_RandomElmPosition(const int zone, const int elm, double *p)
{
  // random numbers
  const double
    rand1 = plas_RandomDouble(),
    rand2 = plas_RandomDouble(),
    rand3 = plas_RandomDouble(),
    rand4 = plas_RandomDouble();

  // get element nodes
  std::vector< unsigned > enu;
  std::vector< int      > en;
  getElmNodes(zone,elm,enu);
  en.resize(enu.size());
  for (size_t i=0; i<enu.size(); ++i)
    en[i] = (int) enu[i];

  switch (getElmType(zone,elm)) {
  case ELM_TRIANGLE: {

    for (int d=0; d<fp.numDim; ++d)
      p[d] =                   plas_getQuantity(COORD_X+d,en[0])
                      + rand1*(plas_getQuantity(COORD_X+d,en[1]) - plas_getQuantity(COORD_X+d,en[0]))
           + (1.-rand1)*rand2*(plas_getQuantity(COORD_X+d,en[2]) - plas_getQuantity(COORD_X+d,en[0]));

  } break;
  case ELM_QUAD: {

    for (int d=0; d<fp.numDim; ++d)
      p[d] = rand3*(                                    plas_getQuantity(COORD_X+d,en[0])
           + rand1*(plas_getQuantity(COORD_X+d,en[1]) - plas_getQuantity(COORD_X+d,en[0])))
           + (1.-rand3) * (                             plas_getQuantity(COORD_X+d,en[2])
           + rand2*(plas_getQuantity(COORD_X+d,en[3]) - plas_getQuantity(COORD_X+d,en[2])));

  } break;
  case ELM_TETRAHEDRON: {

    for (int d=0; d<fp.numDim; ++d)
      p[d] =                                    plas_getQuantity(COORD_X+d,en[0])
                                       + rand1*(plas_getQuantity(COORD_X+d,en[1]) - plas_getQuantity(COORD_X+d,en[0]))
                            + (1.-rand1)*rand2*(plas_getQuantity(COORD_X+d,en[2]) - plas_getQuantity(COORD_X+d,en[0]))
           + (1.-rand1-(1.-rand1)*rand2)*rand3*(plas_getQuantity(COORD_X+d,en[3]) - plas_getQuantity(COORD_X+d,en[0]));

  } break;
  case ELM_HEX: {

    for (int d=0; d<fp.numDim; ++d)
      p[d] = rand4*(rand3*(plas_getQuantity(COORD_X+d,en[0]) + rand1*(plas_getQuantity(COORD_X+d,en[4]) - plas_getQuantity(COORD_X+d,en[0])))
             + (1.-rand3)*(plas_getQuantity(COORD_X+d,en[1]) + rand2*(plas_getQuantity(COORD_X+d,en[5]) - plas_getQuantity(COORD_X+d,en[1]))))
       +(1.-rand4)*(rand3*(plas_getQuantity(COORD_X+d,en[2]) + rand1*(plas_getQuantity(COORD_X+d,en[6]) - plas_getQuantity(COORD_X+d,en[2])))
             + (1.-rand3)*(plas_getQuantity(COORD_X+d,en[3]) + rand2*(plas_getQuantity(COORD_X+d,en[7]) - plas_getQuantity(COORD_X+d,en[3]))));

  } break;
  case ELM_PRISM: {

    for (int d=0; d<fp.numDim; ++d)
      p[d] =                        rand3*plas_getQuantity(COORD_X+d,en[0])
                           + rand3*rand1*(plas_getQuantity(COORD_X+d,en[1]) - plas_getQuantity(COORD_X+d,en[0]))
                + rand3*(1.-rand1)*rand2*(plas_getQuantity(COORD_X+d,en[2]) - plas_getQuantity(COORD_X+d,en[0]))
                             + (1.-rand3)*plas_getQuantity(COORD_X+d,en[3])
                      + (1.-rand3)*rand1*(plas_getQuantity(COORD_X+d,en[4]) - plas_getQuantity(COORD_X+d,en[3]))
           + (1.-rand3)*(1.-rand1)*rand2*(plas_getQuantity(COORD_X+d,en[5]) - plas_getQuantity(COORD_X+d,en[3]));

  } break;
  case ELM_PYRAMID: {

    for (int d=0; d<fp.numDim; ++d)
      p[d] =              rand4*plas_getQuantity(COORD_X+d,en[4])
           + (1.-rand4)*(rand3*(plas_getQuantity(COORD_X+d,en[0])
                       + rand1*(plas_getQuantity(COORD_X+d,en[1])
                              - plas_getQuantity(COORD_X+d,en[0])))
                  + (1.-rand3)*(plas_getQuantity(COORD_X+d,en[2])
                       + rand2*(plas_getQuantity(COORD_X+d,en[3])
                              - plas_getQuantity(COORD_X+d,en[2]))));

  } break;
  default: break;
  }

}


double plas::plas_getElmFaceMiddlePoint(const int iz, const int ie, const int iface, const int dim)
{
  std::vector< int > fn(4,-1);
  plas_getElmFaceNodes(iz,ie,iface,&fn[0],&fn[1],&fn[2],&fn[3]);

  double coord = 0.;
  int ctr = 0;
  for (int i=0; i<4; ++i) {
    if (fn[i]!=-1) {
      coord += plas_getQuantity(COORD_X+dim,fn[i]);
      ++ctr;
    }
  }
  return (ctr? coord/double(ctr) : 0.);
}


void plas::plas_getElmFaceNodes(const int iz, const int ie, const int iface, int *n1, int *n2, int *n3, int *n4)
{
  // get element nodes
  std::vector< unsigned > enu;
  std::vector< int      > en;
  getElmNodes(iz,ie,enu);
  en.resize(enu.size());
  for (size_t i=0; i<enu.size(); ++i)
    en[i] = (int) enu[i];

  *n1 = *n2 = *n3 = *n4 = -1;
  switch (getElmType(iz,ie)) {
  case ELM_TRIANGLE:
    switch (iface) {
    case 0: { *n1=en[0]; *n2=en[1]; } break;
    case 1: { *n1=en[1]; *n2=en[2]; } break;
    case 2: { *n1=en[2]; *n2=en[0]; } break;
    } break;
  case ELM_TETRAHEDRON:
    switch (iface) {
    case 0: { *n1=en[1]; *n2=en[0]; *n3=en[2]; } break;
    case 1: { *n1=en[0]; *n2=en[1]; *n3=en[3]; } break;
    case 2: { *n1=en[1]; *n2=en[2]; *n3=en[3]; } break;
    case 3: { *n1=en[2]; *n2=en[0]; *n3=en[3]; } break;
    } break;
  case ELM_QUAD:
    switch (iface) {
    case 0: { *n1=en[0]; *n2=en[1]; } break;
    case 1: { *n1=en[1]; *n2=en[2]; } break;
    case 2: { *n1=en[2]; *n2=en[3]; } break;
    case 3: { *n1=en[3]; *n2=en[0]; } break;
    } break;
  case ELM_HEX:
    switch (iface) {
    case 0: { *n1=en[0]; *n2=en[1]; *n3=en[5]; *n4=en[4]; } break;
    case 1: { *n1=en[1]; *n2=en[3]; *n3=en[7]; *n4=en[5]; } break;
    case 2: { *n1=en[3]; *n2=en[2]; *n3=en[6]; *n4=en[7]; } break;
    case 3: { *n1=en[2]; *n2=en[0]; *n3=en[4]; *n4=en[6]; } break;
    case 4: { *n1=en[1]; *n2=en[0]; *n3=en[2]; *n4=en[3]; } break;
    case 5: { *n1=en[4]; *n2=en[5]; *n3=en[7]; *n4=en[6]; } break;
    } break;
  case ELM_PRISM:
    switch (iface) {
    case 0: { *n1=en[0]; *n2=en[1]; *n3=en[4]; *n4=en[3]; } break;
    case 1: { *n1=en[1]; *n2=en[2]; *n3=en[5]; *n4=en[4]; } break;
    case 2: { *n1=en[2]; *n2=en[0]; *n3=en[3]; *n4=en[5]; } break;
    case 3: { *n1=en[0]; *n2=en[2]; *n3=en[1]; } break;
    case 4: { *n1=en[3]; *n2=en[4]; *n3=en[5]; } break;
    } break;
  case ELM_PYRAMID:
    switch (iface) {
    case 0: { *n1=en[0]; *n2=en[2]; *n3=en[3]; *n4=en[1]; } break;
    case 1: { *n1=en[0]; *n2=en[1]; *n3=en[4]; } break;
    case 2: { *n1=en[1]; *n2=en[3]; *n3=en[4]; } break;
    case 3: { *n1=en[3]; *n2=en[2]; *n3=en[4]; } break;
    case 4: { *n1=en[2]; *n2=en[0]; *n3=en[4]; } break;
    } break;
  }
}


int plas::plas_getElmNNodes(const int iz, const int ie)
{
  const int t(getElmType(iz,ie));
  return (t==ELM_TRIANGLE?    3 :
         (t==ELM_TETRAHEDRON? 4 :
         (t==ELM_PRISM?       6 :
         (t==ELM_QUAD?        4 :
         (t==ELM_HEX?         8 :
         (t==ELM_PYRAMID?     5 : 0 ))))));
}

