
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


double plas::plas_CalcElmSize(const unsigned iz, const unsigned ie)
{
  plas_elmtype_t et(getElmType(iz,ie));
  std::vector< int > n(plas_getElmNNodes(et),-1);
  getElmNodes(iz,ie,&n[0]);

  switch (et) {
  case ELM_TRIANGLE:    return plas_CalcSizeTriangle(n[0],n[1],n[2]);
  case ELM_TETRAHEDRON: return plas_CalcSizeTetrahedron(n[0],n[1],n[2],n[3]);
  case ELM_WEDGE:       return plas_CalcSizeTetrahedron(n[0],n[1],n[2],n[5])
                             + plas_CalcSizeTetrahedron(n[0],n[1],n[5],n[4])
                             + plas_CalcSizeTetrahedron(n[0],n[4],n[5],n[3]);
  case ELM_QUAD:        return plas_CalcSizeTriangle(n[0],n[1],n[2])
                             + plas_CalcSizeTriangle(n[0],n[2],n[3]);
  case ELM_BRICK:       return plas_CalcSizeTetrahedron(n[0],n[1],n[3],n[5])
                             + plas_CalcSizeTetrahedron(n[0],n[3],n[6],n[5])
                             + plas_CalcSizeTetrahedron(n[0],n[3],n[2],n[6])
                             + plas_CalcSizeTetrahedron(n[0],n[5],n[6],n[4])
                             + plas_CalcSizeTetrahedron(n[3],n[6],n[5],n[7]);
  case ELM_PYRAMID:     return plas_CalcSizeTetrahedron(n[0],n[1],n[3],n[4])
                             + plas_CalcSizeTetrahedron(n[0],n[3],n[2],n[4]);
  case ELM_EDGE:
  case ELM_UNDEFINED:
  case ALL_ELEMENTS:
  default: break;
  }
  return 0.;
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
  std::vector< int > en(plas_getElmNNodes(getElmType(zone,elm)),-1);
  getElmNodes(zone,elm,&en[0]);

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
  case ELM_BRICK: {

    for (int d=0; d<fp.numDim; ++d)
      p[d] = rand4*(rand3*(plas_getQuantity(COORD_X+d,en[0]) + rand1*(plas_getQuantity(COORD_X+d,en[4]) - plas_getQuantity(COORD_X+d,en[0])))
             + (1.-rand3)*(plas_getQuantity(COORD_X+d,en[1]) + rand2*(plas_getQuantity(COORD_X+d,en[5]) - plas_getQuantity(COORD_X+d,en[1]))))
       +(1.-rand4)*(rand3*(plas_getQuantity(COORD_X+d,en[2]) + rand1*(plas_getQuantity(COORD_X+d,en[6]) - plas_getQuantity(COORD_X+d,en[2])))
             + (1.-rand3)*(plas_getQuantity(COORD_X+d,en[3]) + rand2*(plas_getQuantity(COORD_X+d,en[7]) - plas_getQuantity(COORD_X+d,en[3]))));

  } break;
  case ELM_WEDGE: {

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


void plas::plas_getElmFaceMiddlePoint(const int iz, const int ie, const int iface, double *fmp)
{
  plas_elmtype_t ft = plas_getElmFaceType(getElmType(iz,ie),iface);
  std::vector< int > fn(plas_getElmNNodes(ft),-1);
  plas_getElmFaceNodes(iz,ie,iface,&fn[0]);

  for (int d=0; d<fp.numDim; ++d) {
    fmp[d] = 0.;
    for (int i=0; i<plas_getElmNNodes(ft); ++i)
      fmp[d] += plas_getQuantity(COORD_X+d,fn[i]);
    fmp[d] /= (double) plas_getElmNNodes(ft);
  }
}


plas_elmtype_t plas::plas_getElmFaceType(const plas_elmtype_t et, const int iface)
{
  return (et==ELM_TRIANGLE?    ELM_EDGE     :
         (et==ELM_TETRAHEDRON? ELM_TRIANGLE :
         (et==ELM_QUAD?        ELM_EDGE     :
         (et==ELM_BRICK?       ELM_QUAD     :
         (et==ELM_WEDGE?       (iface<3? ELM_QUAD : ELM_TRIANGLE) :
         (et==ELM_PYRAMID?     (iface<1? ELM_QUAD : ELM_TRIANGLE) :
                               ELM_UNDEFINED ))))));
}


void plas::plas_getElmFaceNodes(const int iz, const int ie, const int iface, int *fnodes)
{
  // get element nodes
  const plas_elmtype_t et = getElmType(iz,ie);
  std::vector< int > en(plas_getElmNNodes(et),-1);
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


int plas::plas_getElmNFaces(const plas_elmtype_t et)
{
  return (et==ELM_TRIANGLE?    3 :
         (et==ELM_TETRAHEDRON? 4 :
         (et==ELM_WEDGE?       5 :
         (et==ELM_QUAD?        4 :
         (et==ELM_BRICK?       6 :
         (et==ELM_PYRAMID?     5 : 0 ))))));
}


int plas::plas_getElmNNodes(const plas_elmtype_t et)
{
  return (et==ELM_TRIANGLE?    3 :
         (et==ELM_TETRAHEDRON? 4 :
         (et==ELM_WEDGE?       6 :
         (et==ELM_QUAD?        4 :
         (et==ELM_BRICK?       8 :
         (et==ELM_PYRAMID?     5 : 0 ))))));
}

