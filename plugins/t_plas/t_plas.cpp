
#include <numeric>
#include "ext/xmlParser.h"
#include "mfactory.h"
#include "t_plas.h"

using namespace m;


Register< mtransform,t_plas > mt_plas( 7,
  "-tplas", "[str] filename or string with xml formatted as:",
  "", "<plas",
  "", " iterations=\"[int]\" (default 1)",
  "", " dt=\"[real]\" (default 1.)",
  "", " material.continuum=\"(air|water|nitrogen)\" (default water)",
  "", " <wall zone=\"[str]\"/> (walls, 1 or more)",
  "", "</plas>" );


namespace t_plas_aux {


  unsigned getzoneidx(const std::string& n, const m::mmesh& m)
  {
    for (unsigned r=0; r!=m.z(); ++r)
      if (m.vz[r].n==n)
        return r;
    std::cerr << "error: zone \"" << n << "\" not present!" << std::endl;
    throw 42;
    return 0;
  }


  std::vector< m::mzone >::const_iterator getzoneit(const std::string& n, const m::mmesh& m)
  {
    return m.vz.begin() + getzoneidx(n,m);
  }


}


void t_plas::transform(GetPot& o, mmesh& m)
{
  using std::cout;
  using std::endl;


  cout << "info: setup plas xml..." << endl;
  const std::string o_xml = o.get(o.inc_cursor(),"");
  XMLNode x = ((o_xml.size()? o_xml[0]:'?')=='<')? XMLNode::parseString(o_xml.c_str(),"plas")
                                                 : XMLNode::openFileHelper(o_xml.c_str(),"plas");

  dparam.numIter =  x.getAttribute< int    >("iterations",1 );
  dparam.dt      =  x.getAttribute< double >("dt",        1.);
  const std::string mat =  x.getAttribute< std::string >("continuum.material","water");

  dparam.numIter  = dparam.numIter<0? 0      : dparam.numIter;
  dparam.dt       = dparam.dt<=0.?    1.e-12 : dparam.dt;
  dparam.material = (mat=="air"?      AIR      :
                    (mat=="water"?    WATER    :
                    (mat=="nitrogen"? NITROGEN : WATER    )));
  cout << "info: setup plas xml." << endl;


  screenOutput("converting m::mmesh to dmesh...");
  M = &m;
  plasdriver_ReadGambitNeutralFile(x);
  plasdriver_CalcElmsAroundNode();
  plasdriver_CalcElementNeighbours();
  plasdriver_CalcElementNormals();
  plasdriver_CalcElementVolumes();
  plasdriver_CalcNodalVolumes();


  screenOutput("initializing flow field...");
  dflow.p = new double [M->n()];
  dflow.T = new double [M->n()];
  dflow.u = new double*[M->n()];
  for (unsigned i=0; i<M->n(); ++i)
    dflow.u[i] = new double[M->d()];
  plasdriver_InitFlowField(dparam.material);


  screenOutput("initializing PLaS...");
  plas::initialize(x);
  screenOutput("initializing PLaS.");


  screenOutput("perform PLaS iterations...");
  for (dparam.iter=1; dparam.iter<=dparam.numIter; ++dparam.iter) {
    screenOutput("iterate...");
    plas::run();
  }
  screenOutput("perform PLaS iterations.");


  screenOutput("terminating PLaS...");
  ///plasdriver_FreeGambitMemory();
  for (unsigned i=0; i<M->n(); ++i)
    delete[] dflow.u[i];
  delete[] dflow.p;
  delete[] dflow.u;
  delete[] dflow.T;
  screenOutput("terminating PLaS.");
}


// This file contains routines to compute the geometry of the
// mesh used for the steady-state flow solution.
// This function computes the area of a triangle.
double t_plas::plasdriver_CalcAreaTriangle(double c[3][2])
{
  return c[0][0]*(c[1][1]-c[2][1]) + c[1][0]*(c[2][1]-c[0][1]) + c[2][0]*(c[0][1]-c[1][1]);
}


// This file contains routines to compute the geometry of the
// mesh used for the steady-state flow solution.
// This function computes element neighbours.
void t_plas::plasdriver_CalcElementNeighbours()
{
  screenOutput("performing geometry calculations: element neighbours...");

  int neighbourFound,jnod,knod,lnod,ielm,jelm,kelm,lelm,ifac,faceNodes[4];

  dmesh.elmNeighbs.resize(dmesh.numElm);
  for(ielm=0; ielm<dmesh.numElm; ielm++){
    dmesh.elmNeighbs[ielm].assign(dmesh.numElmFaces[ielm],0);
    for(ifac=0; ifac<dmesh.numElmFaces[ielm]; ifac++){
      plasdriver_GetFaceNodes(ielm,ifac,faceNodes);
      neighbourFound = 0;
      for(jnod=0; jnod<dmesh.numNodElms[faceNodes[0]]; jnod++){
        jelm = dmesh.nodElms[faceNodes[0]][jnod];
        for(knod=0; knod<dmesh.numNodElms[faceNodes[1]]; knod++){
          kelm = dmesh.nodElms[faceNodes[1]][knod];
          if(M->d()==2){
            if(jelm==kelm && jelm!=ielm){neighbourFound = 1;}
          }else if(M->d()==3){
            for(lnod=0; lnod<dmesh.numNodElms[faceNodes[2]]; lnod++){
              lelm = dmesh.nodElms[faceNodes[2]][lnod];
              if(jelm==kelm && jelm==lelm && jelm!=ielm){neighbourFound = 1;}
              if(neighbourFound==1){break;}
            }
          }
          if(neighbourFound==1){break;}
        }
        if(neighbourFound==1){break;}
      }
      if(neighbourFound==1){
        dmesh.elmNeighbs[ielm][ifac] = jelm;
      } else{
        dmesh.elmNeighbs[ielm][ifac] = -1;
      }
    }
  }
}


// This file contains routines to compute the geometry of the
// mesh used for the steady-state flow solution.
// This function computes element normals.
void t_plas::plasdriver_CalcElementNormals()
{
  screenOutput("performing geometry calculations: element normals...");

  int ielm,ifac,fnodes[4];

  for(ielm=0; ielm<dmesh.numElm; ielm++){
    for(ifac=0; ifac<dmesh.numElmFaces[ielm]; ifac++){
      plasdriver_GetFaceNodes(ielm,ifac,fnodes);

      if(dmesh.elmTypes[ielm]==ELM_SIMPLEX && M->d()==2){

        dmesh.elmNorms[ielm][ifac][0] = M->vv[1][fnodes[0]] - M->vv[1][fnodes[1]];
        dmesh.elmNorms[ielm][ifac][1] = M->vv[0][fnodes[1]] - M->vv[0][fnodes[0]];

      } else if((dmesh.elmTypes[ielm]==ELM_SIMPLEX && M->d()==3)
                || (dmesh.elmTypes[ielm]==ELM_PRISM && ifac>2)
                || (dmesh.elmTypes[ielm]==ELM_PYRAMID && ifac>0)){

        dmesh.elmNorms[ielm][ifac][0] =
          0.5*((M->vv[1][fnodes[2]] - M->vv[1][fnodes[0]])
              *(M->vv[2][fnodes[1]] - M->vv[2][fnodes[0]])
              -(M->vv[1][fnodes[1]] - M->vv[1][fnodes[0]])
              *(M->vv[2][fnodes[2]] - M->vv[2][fnodes[0]]));
        dmesh.elmNorms[ielm][ifac][1] =
          0.5*((M->vv[2][fnodes[2]] - M->vv[2][fnodes[0]])
              *(M->vv[0][fnodes[1]] - M->vv[0][fnodes[0]])
              -(M->vv[2][fnodes[1]] - M->vv[2][fnodes[0]])
              *(M->vv[0][fnodes[2]] - M->vv[0][fnodes[0]]));
        dmesh.elmNorms[ielm][ifac][2] =
          0.5*((M->vv[0][fnodes[2]] - M->vv[0][fnodes[0]])
              *(M->vv[1][fnodes[1]] - M->vv[1][fnodes[0]])
              -(M->vv[0][fnodes[1]] - M->vv[0][fnodes[0]])
              *(M->vv[1][fnodes[2]] - M->vv[1][fnodes[0]]));

      } else if(dmesh.elmTypes[ielm]==ELM_QUAD){

        dmesh.elmNorms[ielm][ifac][0] = 2.0*(M->vv[1][fnodes[0]] - M->vv[1][fnodes[1]]);
        dmesh.elmNorms[ielm][ifac][1] = 2.0*(M->vv[0][fnodes[1]] - M->vv[0][fnodes[0]]);

      } else if(dmesh.elmTypes[ielm]==ELM_HEX
                || (dmesh.elmTypes[ielm]==ELM_PRISM && ifac<=2)
                || (dmesh.elmTypes[ielm]==ELM_PYRAMID && ifac==0)){

        dmesh.elmNorms[ielm][ifac][0] =
            ((M->vv[1][fnodes[3]] - M->vv[1][fnodes[0]])
            *(M->vv[2][fnodes[1]] - M->vv[2][fnodes[0]])
            -(M->vv[1][fnodes[1]] - M->vv[1][fnodes[0]])
            *(M->vv[2][fnodes[3]] - M->vv[2][fnodes[0]]));
        dmesh.elmNorms[ielm][ifac][1] =
            ((M->vv[2][fnodes[3]] - M->vv[2][fnodes[0]])
            *(M->vv[0][fnodes[1]] - M->vv[0][fnodes[0]])
            -(M->vv[2][fnodes[1]] - M->vv[2][fnodes[0]])
            *(M->vv[0][fnodes[3]] - M->vv[0][fnodes[0]]));
        dmesh.elmNorms[ielm][ifac][2] =
            ((M->vv[0][fnodes[3]] - M->vv[0][fnodes[0]])
            *(M->vv[1][fnodes[1]] - M->vv[1][fnodes[0]])
            -(M->vv[0][fnodes[1]] - M->vv[0][fnodes[0]])
            *(M->vv[1][fnodes[3]] - M->vv[1][fnodes[0]]));
      }
    }
  }
}


// This file contains routines to compute the geometry of the
// mesh used for the steady-state flow solution.
// This function computes element volumes.
void t_plas::plasdriver_CalcElementVolumes()
{
  screenOutput("performing geometry calculations: element volumes...");

  int ielm;
  double c2[3][2],c3[4][3];

  dmesh.elmVolumes.assign(dmesh.numElm,0.);

  for(ielm=0; ielm<dmesh.numElm; ielm++){

    if(dmesh.elmTypes[ielm]==ELM_SIMPLEX && M->d()==2){
      c2[0][0] = M->vv[0][ M->vz[0].e2n[ielm].n[0] ];
      c2[0][1] = M->vv[1][ M->vz[0].e2n[ielm].n[0] ];
      c2[1][0] = M->vv[0][ M->vz[0].e2n[ielm].n[1] ];
      c2[1][1] = M->vv[1][ M->vz[0].e2n[ielm].n[1] ];
      c2[2][0] = M->vv[0][ M->vz[0].e2n[ielm].n[2] ];
      c2[2][1] = M->vv[1][ M->vz[0].e2n[ielm].n[2] ];
      dmesh.elmVolumes[ielm] = plasdriver_CalcAreaTriangle(c2);
    } else if(dmesh.elmTypes[ielm]==ELM_SIMPLEX && M->d()==3){
      c3[0][0] = M->vv[0][ M->vz[0].e2n[ielm].n[0] ];
      c3[0][1] = M->vv[1][ M->vz[0].e2n[ielm].n[0] ];
      c3[0][2] = M->vv[2][ M->vz[0].e2n[ielm].n[0] ];
      c3[1][0] = M->vv[0][ M->vz[0].e2n[ielm].n[1] ];
      c3[1][1] = M->vv[1][ M->vz[0].e2n[ielm].n[1] ];
      c3[1][2] = M->vv[2][ M->vz[0].e2n[ielm].n[1] ];
      c3[2][0] = M->vv[0][ M->vz[0].e2n[ielm].n[2] ];
      c3[2][1] = M->vv[1][ M->vz[0].e2n[ielm].n[2] ];
      c3[2][2] = M->vv[2][ M->vz[0].e2n[ielm].n[2] ];
      c3[3][0] = M->vv[0][ M->vz[0].e2n[ielm].n[3] ];
      c3[3][1] = M->vv[1][ M->vz[0].e2n[ielm].n[3] ];
      c3[3][2] = M->vv[2][ M->vz[0].e2n[ielm].n[3] ];
      dmesh.elmVolumes[ielm] = plasdriver_CalcVolumeTetra(c3);
    } else if(dmesh.elmTypes[ielm]==ELM_QUAD){
      c2[0][0] = M->vv[0][ M->vz[0].e2n[ielm].n[0] ];
      c2[0][1] = M->vv[1][ M->vz[0].e2n[ielm].n[0] ];
      c2[1][0] = M->vv[0][ M->vz[0].e2n[ielm].n[1] ];
      c2[1][1] = M->vv[1][ M->vz[0].e2n[ielm].n[1] ];
      c2[2][0] = M->vv[0][ M->vz[0].e2n[ielm].n[2] ];
      c2[2][1] = M->vv[1][ M->vz[0].e2n[ielm].n[2] ];
      dmesh.elmVolumes[ielm] = plasdriver_CalcAreaTriangle(c2);
      c2[0][0] = M->vv[0][ M->vz[0].e2n[ielm].n[0] ];
      c2[0][1] = M->vv[1][ M->vz[0].e2n[ielm].n[0] ];
      c2[1][0] = M->vv[0][ M->vz[0].e2n[ielm].n[2] ];
      c2[1][1] = M->vv[1][ M->vz[0].e2n[ielm].n[2] ];
      c2[2][0] = M->vv[0][ M->vz[0].e2n[ielm].n[3] ];
      c2[2][1] = M->vv[1][ M->vz[0].e2n[ielm].n[3] ];
      dmesh.elmVolumes[ielm] += plasdriver_CalcAreaTriangle(c2);
    } else if(dmesh.elmTypes[ielm]==ELM_HEX){
      c3[0][0] = M->vv[0][ M->vz[0].e2n[ielm].n[0] ];
      c3[0][1] = M->vv[1][ M->vz[0].e2n[ielm].n[0] ];
      c3[0][2] = M->vv[2][ M->vz[0].e2n[ielm].n[0] ];
      c3[1][0] = M->vv[0][ M->vz[0].e2n[ielm].n[1] ];
      c3[1][1] = M->vv[1][ M->vz[0].e2n[ielm].n[1] ];
      c3[1][2] = M->vv[2][ M->vz[0].e2n[ielm].n[1] ];
      c3[2][0] = M->vv[0][ M->vz[0].e2n[ielm].n[3] ];
      c3[2][1] = M->vv[1][ M->vz[0].e2n[ielm].n[3] ];
      c3[2][2] = M->vv[2][ M->vz[0].e2n[ielm].n[3] ];
      c3[3][0] = M->vv[0][ M->vz[0].e2n[ielm].n[5] ];
      c3[3][1] = M->vv[1][ M->vz[0].e2n[ielm].n[5] ];
      c3[3][2] = M->vv[2][ M->vz[0].e2n[ielm].n[5] ];
      dmesh.elmVolumes[ielm] = plasdriver_CalcVolumeTetra(c3);
      c3[0][0] = M->vv[0][ M->vz[0].e2n[ielm].n[0] ];
      c3[0][1] = M->vv[1][ M->vz[0].e2n[ielm].n[0] ];
      c3[0][2] = M->vv[2][ M->vz[0].e2n[ielm].n[0] ];
      c3[1][0] = M->vv[0][ M->vz[0].e2n[ielm].n[3] ];
      c3[1][1] = M->vv[1][ M->vz[0].e2n[ielm].n[3] ];
      c3[1][2] = M->vv[2][ M->vz[0].e2n[ielm].n[3] ];
      c3[2][0] = M->vv[0][ M->vz[0].e2n[ielm].n[6] ];
      c3[2][1] = M->vv[1][ M->vz[0].e2n[ielm].n[6] ];
      c3[2][2] = M->vv[2][ M->vz[0].e2n[ielm].n[6] ];
      c3[3][0] = M->vv[0][ M->vz[0].e2n[ielm].n[5] ];
      c3[3][1] = M->vv[1][ M->vz[0].e2n[ielm].n[5] ];
      c3[3][2] = M->vv[2][ M->vz[0].e2n[ielm].n[5] ];
      dmesh.elmVolumes[ielm] += plasdriver_CalcVolumeTetra(c3);
      c3[0][0] = M->vv[0][ M->vz[0].e2n[ielm].n[0] ];
      c3[0][1] = M->vv[1][ M->vz[0].e2n[ielm].n[0] ];
      c3[0][2] = M->vv[2][ M->vz[0].e2n[ielm].n[0] ];
      c3[1][0] = M->vv[0][ M->vz[0].e2n[ielm].n[3] ];
      c3[1][1] = M->vv[1][ M->vz[0].e2n[ielm].n[3] ];
      c3[1][2] = M->vv[2][ M->vz[0].e2n[ielm].n[3] ];
      c3[2][0] = M->vv[0][ M->vz[0].e2n[ielm].n[2] ];
      c3[2][1] = M->vv[1][ M->vz[0].e2n[ielm].n[2] ];
      c3[2][2] = M->vv[2][ M->vz[0].e2n[ielm].n[2] ];
      c3[3][0] = M->vv[0][ M->vz[0].e2n[ielm].n[6] ];
      c3[3][1] = M->vv[1][ M->vz[0].e2n[ielm].n[6] ];
      c3[3][2] = M->vv[2][ M->vz[0].e2n[ielm].n[6] ];
      dmesh.elmVolumes[ielm] += plasdriver_CalcVolumeTetra(c3);
      c3[0][0] = M->vv[0][ M->vz[0].e2n[ielm].n[0] ];
      c3[0][1] = M->vv[1][ M->vz[0].e2n[ielm].n[0] ];
      c3[0][2] = M->vv[2][ M->vz[0].e2n[ielm].n[0] ];
      c3[1][0] = M->vv[0][ M->vz[0].e2n[ielm].n[5] ];
      c3[1][1] = M->vv[1][ M->vz[0].e2n[ielm].n[5] ];
      c3[1][2] = M->vv[2][ M->vz[0].e2n[ielm].n[5] ];
      c3[2][0] = M->vv[0][ M->vz[0].e2n[ielm].n[6] ];
      c3[2][1] = M->vv[1][ M->vz[0].e2n[ielm].n[6] ];
      c3[2][2] = M->vv[2][ M->vz[0].e2n[ielm].n[6] ];
      c3[3][0] = M->vv[0][ M->vz[0].e2n[ielm].n[4] ];
      c3[3][1] = M->vv[1][ M->vz[0].e2n[ielm].n[4] ];
      c3[3][2] = M->vv[2][ M->vz[0].e2n[ielm].n[4] ];
      dmesh.elmVolumes[ielm] += plasdriver_CalcVolumeTetra(c3);
      c3[0][0] = M->vv[0][ M->vz[0].e2n[ielm].n[3] ];
      c3[0][1] = M->vv[1][ M->vz[0].e2n[ielm].n[3] ];
      c3[0][2] = M->vv[2][ M->vz[0].e2n[ielm].n[3] ];
      c3[1][0] = M->vv[0][ M->vz[0].e2n[ielm].n[6] ];
      c3[1][1] = M->vv[1][ M->vz[0].e2n[ielm].n[6] ];
      c3[1][2] = M->vv[2][ M->vz[0].e2n[ielm].n[6] ];
      c3[2][0] = M->vv[0][ M->vz[0].e2n[ielm].n[5] ];
      c3[2][1] = M->vv[1][ M->vz[0].e2n[ielm].n[5] ];
      c3[2][2] = M->vv[2][ M->vz[0].e2n[ielm].n[5] ];
      c3[3][0] = M->vv[0][ M->vz[0].e2n[ielm].n[7] ];
      c3[3][1] = M->vv[1][ M->vz[0].e2n[ielm].n[7] ];
      c3[3][2] = M->vv[2][ M->vz[0].e2n[ielm].n[7] ];
      dmesh.elmVolumes[ielm] += plasdriver_CalcVolumeTetra(c3);
    } else if(dmesh.elmTypes[ielm]==ELM_PRISM){
      c3[0][0] = M->vv[0][ M->vz[0].e2n[ielm].n[0] ];
      c3[0][1] = M->vv[1][ M->vz[0].e2n[ielm].n[0] ];
      c3[0][2] = M->vv[2][ M->vz[0].e2n[ielm].n[0] ];
      c3[1][0] = M->vv[0][ M->vz[0].e2n[ielm].n[1] ];
      c3[1][1] = M->vv[1][ M->vz[0].e2n[ielm].n[1] ];
      c3[1][2] = M->vv[2][ M->vz[0].e2n[ielm].n[1] ];
      c3[2][0] = M->vv[0][ M->vz[0].e2n[ielm].n[2] ];
      c3[2][1] = M->vv[1][ M->vz[0].e2n[ielm].n[2] ];
      c3[2][2] = M->vv[2][ M->vz[0].e2n[ielm].n[2] ];
      c3[3][0] = M->vv[0][ M->vz[0].e2n[ielm].n[5] ];
      c3[3][1] = M->vv[1][ M->vz[0].e2n[ielm].n[5] ];
      c3[3][2] = M->vv[2][ M->vz[0].e2n[ielm].n[5] ];
      dmesh.elmVolumes[ielm] = plasdriver_CalcVolumeTetra(c3);
      c3[0][0] = M->vv[0][ M->vz[0].e2n[ielm].n[0] ];
      c3[0][1] = M->vv[1][ M->vz[0].e2n[ielm].n[0] ];
      c3[0][2] = M->vv[2][ M->vz[0].e2n[ielm].n[0] ];
      c3[1][0] = M->vv[0][ M->vz[0].e2n[ielm].n[1] ];
      c3[1][1] = M->vv[1][ M->vz[0].e2n[ielm].n[1] ];
      c3[1][2] = M->vv[2][ M->vz[0].e2n[ielm].n[1] ];
      c3[2][0] = M->vv[0][ M->vz[0].e2n[ielm].n[5] ];
      c3[2][1] = M->vv[1][ M->vz[0].e2n[ielm].n[5] ];
      c3[2][2] = M->vv[2][ M->vz[0].e2n[ielm].n[5] ];
      c3[3][0] = M->vv[0][ M->vz[0].e2n[ielm].n[4] ];
      c3[3][1] = M->vv[1][ M->vz[0].e2n[ielm].n[4] ];
      c3[3][2] = M->vv[2][ M->vz[0].e2n[ielm].n[4] ];
      dmesh.elmVolumes[ielm] += plasdriver_CalcVolumeTetra(c3);
      c3[0][0] = M->vv[0][ M->vz[0].e2n[ielm].n[0] ];
      c3[0][1] = M->vv[1][ M->vz[0].e2n[ielm].n[0] ];
      c3[0][2] = M->vv[2][ M->vz[0].e2n[ielm].n[0] ];
      c3[1][0] = M->vv[0][ M->vz[0].e2n[ielm].n[4] ];
      c3[1][1] = M->vv[1][ M->vz[0].e2n[ielm].n[4] ];
      c3[1][2] = M->vv[2][ M->vz[0].e2n[ielm].n[4] ];
      c3[2][0] = M->vv[0][ M->vz[0].e2n[ielm].n[5] ];
      c3[2][1] = M->vv[1][ M->vz[0].e2n[ielm].n[5] ];
      c3[2][2] = M->vv[2][ M->vz[0].e2n[ielm].n[5] ];
      c3[3][0] = M->vv[0][ M->vz[0].e2n[ielm].n[3] ];
      c3[3][1] = M->vv[1][ M->vz[0].e2n[ielm].n[3] ];
      c3[3][2] = M->vv[2][ M->vz[0].e2n[ielm].n[3] ];
      dmesh.elmVolumes[ielm] += plasdriver_CalcVolumeTetra(c3);
    } else if(dmesh.elmTypes[ielm]==ELM_PYRAMID){
      c3[0][0] = M->vv[0][ M->vz[0].e2n[ielm].n[0] ];
      c3[0][1] = M->vv[1][ M->vz[0].e2n[ielm].n[0] ];
      c3[0][2] = M->vv[2][ M->vz[0].e2n[ielm].n[0] ];
      c3[1][0] = M->vv[0][ M->vz[0].e2n[ielm].n[1] ];
      c3[1][1] = M->vv[1][ M->vz[0].e2n[ielm].n[1] ];
      c3[1][2] = M->vv[2][ M->vz[0].e2n[ielm].n[1] ];
      c3[2][0] = M->vv[0][ M->vz[0].e2n[ielm].n[3] ];
      c3[2][1] = M->vv[1][ M->vz[0].e2n[ielm].n[3] ];
      c3[2][2] = M->vv[2][ M->vz[0].e2n[ielm].n[3] ];
      c3[3][0] = M->vv[0][ M->vz[0].e2n[ielm].n[4] ];
      c3[3][1] = M->vv[1][ M->vz[0].e2n[ielm].n[4] ];
      c3[3][2] = M->vv[2][ M->vz[0].e2n[ielm].n[4] ];
      dmesh.elmVolumes[ielm] = plasdriver_CalcVolumeTetra(c3);
      c3[0][0] = M->vv[0][ M->vz[0].e2n[ielm].n[0] ];
      c3[0][1] = M->vv[1][ M->vz[0].e2n[ielm].n[0] ];
      c3[0][2] = M->vv[2][ M->vz[0].e2n[ielm].n[0] ];
      c3[1][0] = M->vv[0][ M->vz[0].e2n[ielm].n[3] ];
      c3[1][1] = M->vv[1][ M->vz[0].e2n[ielm].n[3] ];
      c3[1][2] = M->vv[2][ M->vz[0].e2n[ielm].n[3] ];
      c3[2][0] = M->vv[0][ M->vz[0].e2n[ielm].n[2] ];
      c3[2][1] = M->vv[1][ M->vz[0].e2n[ielm].n[2] ];
      c3[2][2] = M->vv[2][ M->vz[0].e2n[ielm].n[2] ];
      c3[3][0] = M->vv[0][ M->vz[0].e2n[ielm].n[4] ];
      c3[3][1] = M->vv[1][ M->vz[0].e2n[ielm].n[4] ];
      c3[3][2] = M->vv[2][ M->vz[0].e2n[ielm].n[4] ];
      dmesh.elmVolumes[ielm] += plasdriver_CalcVolumeTetra(c3);
    }
  }

  dmesh.minElmVolume = dmesh.elmVolumes[0];
  dmesh.maxElmVolume = dmesh.elmVolumes[0];
  for (int i=0; i<dmesh.numElm; ++i){
    if (dmesh.elmVolumes[i]<dmesh.minElmVolume) dmesh.minElmVolume = dmesh.elmVolumes[i];
    if (dmesh.elmVolumes[i]>dmesh.maxElmVolume) dmesh.maxElmVolume = dmesh.elmVolumes[i];
  }
}


// This file contains routines to compute the geometry of the
// mesh used for the steady-state flow solution.
// This function computes the elements around a node.
void t_plas::plasdriver_CalcElmsAroundNode()
{
  screenOutput("performing geometry calculations: elements around nodes...");

  dmesh.numNodElms.assign(M->n(),0);
  dmesh.nodElms   .assign(M->n(),std::vector< int >(50,0));  // FIXME

  for (int ielm=0; ielm<dmesh.numElm; ++ielm) {
    for (int inod=0; inod<dmesh.numElmNodes[ielm]; ++inod) {
      dmesh.nodElms[M->vz[0].e2n[ielm].n[inod]][dmesh.numNodElms[M->vz[0].e2n[ielm].n[inod]]] = ielm;
      dmesh.numNodElms[M->vz[0].e2n[ielm].n[inod]]++;
    }
  }
}


// This file contains routines to compute the geometry of the
// mesh used for the steady-state flow solution.
// This function computes the nodal dual cell volumes.
void t_plas::plasdriver_CalcNodalVolumes()
{
  screenOutput("performing geometry calculations: nodal volumes...");

  dmesh.nodVolumes.assign(M->n(),0.);
  for (unsigned n=0; n<M->n(); ++n)
    for (int e=0; e<dmesh.numNodElms[n]; ++e)
      dmesh.nodVolumes[n] += dmesh.elmVolumes[dmesh.nodElms[n][e]]/dmesh.numNodElms[n];
}


// This file contains routines to compute the geometry of the
// mesh used for the steady-state flow solution.
// This function computes the volume of a tetrahedron.
double t_plas::plasdriver_CalcVolumeTetra(double c[4][3])
{
  double v1[3],v2[3],v3[3],v4[3];

  v1[0] = c[1][0]-c[0][0];
  v1[1] = c[1][1]-c[0][1];
  v1[2] = c[1][2]-c[0][2];
  v2[0] = c[2][0]-c[0][0];
  v2[1] = c[2][1]-c[0][1];
  v2[2] = c[2][2]-c[0][2];

  plas_CalcCrossProduct_3D(v3,v1,v2);

  v4[0] = c[3][0]-c[0][0];
  v4[1] = c[3][1]-c[0][1];
  v4[2] = c[3][2]-c[0][2];

  return (plas_CalcVectScalarProduct(3,v3,v4)/6.0);
}


// This file contains routines to compute the geometry of the
// mesh used for the steady-state flow solution.
// This function gets the nodes of a boundary face.
void t_plas::plasdriver_GetFaceNodes(int elm, int face, int *nodes)
{
  if(dmesh.elmTypes[elm]==ELM_SIMPLEX && M->d()==2){
    if(face==0){
      nodes[0] = M->vz[0].e2n[elm].n[0];
      nodes[1] = M->vz[0].e2n[elm].n[1];
      nodes[2] = -1;
      nodes[3] = -1;
    } else if(face==1){
      nodes[0] = M->vz[0].e2n[elm].n[1];
      nodes[1] = M->vz[0].e2n[elm].n[2];
      nodes[2] = -1;
      nodes[3] = -1;
    } else if(face==2){
      nodes[0] = M->vz[0].e2n[elm].n[2];
      nodes[1] = M->vz[0].e2n[elm].n[0];
      nodes[2] = -1;
      nodes[3] = -1;
    }
  } else if(dmesh.elmTypes[elm]==ELM_SIMPLEX && M->d()==3){
    if(face==0){
      nodes[0] = M->vz[0].e2n[elm].n[1];
      nodes[1] = M->vz[0].e2n[elm].n[0];
      nodes[2] = M->vz[0].e2n[elm].n[2];
      nodes[3] = -1;
    } else if(face==1){
      nodes[0] = M->vz[0].e2n[elm].n[0];
      nodes[1] = M->vz[0].e2n[elm].n[1];
      nodes[2] = M->vz[0].e2n[elm].n[3];
      nodes[3] = -1;
    } else if(face==2){
      nodes[0] = M->vz[0].e2n[elm].n[1];
      nodes[1] = M->vz[0].e2n[elm].n[2];
      nodes[2] = M->vz[0].e2n[elm].n[3];
      nodes[3] = -1;
    } else if(face==3){
      nodes[0] = M->vz[0].e2n[elm].n[2];
      nodes[1] = M->vz[0].e2n[elm].n[0];
      nodes[2] = M->vz[0].e2n[elm].n[3];
      nodes[3] = -1;
    }
  } else if(dmesh.elmTypes[elm]==ELM_QUAD){
    if(face==0){
      nodes[0] = M->vz[0].e2n[elm].n[0];
      nodes[1] = M->vz[0].e2n[elm].n[1];
      nodes[2] = -1;
      nodes[3] = -1;
    } else if(face==1){
      nodes[0] = M->vz[0].e2n[elm].n[1];
      nodes[1] = M->vz[0].e2n[elm].n[2];
      nodes[2] = -1;
      nodes[3] = -1;
    } else if(face==2){
      nodes[0] = M->vz[0].e2n[elm].n[2];
      nodes[1] = M->vz[0].e2n[elm].n[3];
      nodes[2] = -1;
      nodes[3] = -1;
    } else if(face==3){
      nodes[0] = M->vz[0].e2n[elm].n[3];
      nodes[1] = M->vz[0].e2n[elm].n[0];
      nodes[2] = -1;
      nodes[3] = -1;
    }
  } else if(dmesh.elmTypes[elm]==ELM_HEX){
    if(face==0){
      nodes[0] = M->vz[0].e2n[elm].n[0];
      nodes[1] = M->vz[0].e2n[elm].n[1];
      nodes[2] = M->vz[0].e2n[elm].n[5];
      nodes[3] = M->vz[0].e2n[elm].n[4];
    } else if(face==1){
      nodes[0] = M->vz[0].e2n[elm].n[1];
      nodes[1] = M->vz[0].e2n[elm].n[3];
      nodes[2] = M->vz[0].e2n[elm].n[7];
      nodes[3] = M->vz[0].e2n[elm].n[5];
    } else if(face==2){
      nodes[0] = M->vz[0].e2n[elm].n[3];
      nodes[1] = M->vz[0].e2n[elm].n[2];
      nodes[2] = M->vz[0].e2n[elm].n[6];
      nodes[3] = M->vz[0].e2n[elm].n[7];
    } else if(face==3){
      nodes[0] = M->vz[0].e2n[elm].n[2];
      nodes[1] = M->vz[0].e2n[elm].n[0];
      nodes[2] = M->vz[0].e2n[elm].n[4];
      nodes[3] = M->vz[0].e2n[elm].n[6];
    } else if(face==4){
      nodes[0] = M->vz[0].e2n[elm].n[1];
      nodes[1] = M->vz[0].e2n[elm].n[0];
      nodes[2] = M->vz[0].e2n[elm].n[2];
      nodes[3] = M->vz[0].e2n[elm].n[3];
    } else if(face==5){
      nodes[0] = M->vz[0].e2n[elm].n[4];
      nodes[1] = M->vz[0].e2n[elm].n[5];
      nodes[2] = M->vz[0].e2n[elm].n[7];
      nodes[3] = M->vz[0].e2n[elm].n[6];
    }
  } else if(dmesh.elmTypes[elm]==ELM_PRISM){
    if(face==0){
      nodes[0] = M->vz[0].e2n[elm].n[0];
      nodes[1] = M->vz[0].e2n[elm].n[1];
      nodes[2] = M->vz[0].e2n[elm].n[4];
      nodes[3] = M->vz[0].e2n[elm].n[3];
    } else if(face==1){
      nodes[0] = M->vz[0].e2n[elm].n[1];
      nodes[1] = M->vz[0].e2n[elm].n[2];
      nodes[2] = M->vz[0].e2n[elm].n[5];
      nodes[3] = M->vz[0].e2n[elm].n[4];
    } else if(face==2){
      nodes[0] = M->vz[0].e2n[elm].n[2];
      nodes[1] = M->vz[0].e2n[elm].n[0];
      nodes[2] = M->vz[0].e2n[elm].n[3];
      nodes[3] = M->vz[0].e2n[elm].n[5];
    } else if(face==3){
      nodes[0] = M->vz[0].e2n[elm].n[0];
      nodes[1] = M->vz[0].e2n[elm].n[2];
      nodes[2] = M->vz[0].e2n[elm].n[1];
      nodes[3] = -1;
    } else if(face==4){
      nodes[0] = M->vz[0].e2n[elm].n[3];
      nodes[1] = M->vz[0].e2n[elm].n[4];
      nodes[2] = M->vz[0].e2n[elm].n[5];
      nodes[3] = -1;
    }
  } else if(dmesh.elmTypes[elm]==ELM_PYRAMID){
    if(face==0){
      nodes[0] = M->vz[0].e2n[elm].n[0];
      nodes[1] = M->vz[0].e2n[elm].n[2];
      nodes[2] = M->vz[0].e2n[elm].n[3];
      nodes[3] = M->vz[0].e2n[elm].n[1];
    } else if(face==1){
      nodes[0] = M->vz[0].e2n[elm].n[0];
      nodes[1] = M->vz[0].e2n[elm].n[1];
      nodes[2] = M->vz[0].e2n[elm].n[4];
      nodes[3] = -1;
    } else if(face==2){
      nodes[0] = M->vz[0].e2n[elm].n[1];
      nodes[1] = M->vz[0].e2n[elm].n[3];
      nodes[2] = M->vz[0].e2n[elm].n[4];
      nodes[3] = -1;
    } else if(face==3){
      nodes[0] = M->vz[0].e2n[elm].n[3];
      nodes[1] = M->vz[0].e2n[elm].n[2];
      nodes[2] = M->vz[0].e2n[elm].n[4];
      nodes[3] = -1;
    } else if(face==4){
      nodes[0] = M->vz[0].e2n[elm].n[2];
      nodes[1] = M->vz[0].e2n[elm].n[0];
      nodes[2] = M->vz[0].e2n[elm].n[4];
      nodes[3] = -1;
    }
  }
}


// This file contains functionality to define or read a
// steady-state primary flow field.
// This function initializes a steady-state flow field.
void t_plas::plasdriver_InitFlowField(int material)
{
  // Define pressure, velocity and temperature
  double
    p = 101325.,
    u =      1.,
    v =      0.,
    w =      0.,
    T =    373.124;

  if (material==AIR) {

    dflow.rho = 1.225;
    dflow.mu  = 1.7894e-5;
    dflow.cp  = 1.006;
    dflow.k   = 0.0242;

  }
  else if (material==WATER) {

    dflow.rho = 998.2;
    dflow.mu  = 1.003e-3;
    dflow.cp  = 4.182;
    dflow.k   = 0.6;

  }
  else if (material==NITROGEN) {

    dflow.rho = 1.25;
    dflow.mu  = (6.5592e-7*pow(T,0.6081))/(1.0+54.715/T);
    dflow.cp  = (6.50+0.001*T)*4.184/(28.01e-3);
    dflow.k   = 2.5*(6.5592e-7*pow(T,0.6081))/(1.0+54.715/T)*((6.50+0.001*T)*4.184/(28.01e-3)-8.314)/28.01;

  }

  // Impose flow variables
  for (unsigned n=0; n<M->n(); ++n) {
    dflow.p[n]    = p;
    dflow.u[n][0] = u;
    dflow.u[n][1] = v;
    dflow.u[n][2] = w;
    dflow.T[n]    = T;
  }
}


// This file contains routines to compute the geometry of the
// mesh used for the steady-state flow solution.
// This function reads a Gambit mesh.
void t_plas::plasdriver_ReadGambitNeutralFile(const XMLNode& x)
{
  // general information
  m_zinner_nelems.assign(M->z(),0);
  for (unsigned i=0; i<M->z(); ++i)
    m_zinner_nelems[i] = (int)(M->d(i)==M->d()? M->e(i) : 0);
  dmesh.numElm = std::accumulate(m_zinner_nelems.begin(),m_zinner_nelems.end(),0);


  // set volume elements
  dmesh.elmTypes   .assign(dmesh.numElm,0);
  dmesh.numElmNodes.assign(dmesh.numElm,0);
  dmesh.numElmFaces.assign(dmesh.numElm,0);
  dmesh.elmNorms   .resize(dmesh.numElm);
  for (unsigned iz=0, ie=0; iz<M->z(); ++iz) {
    if (M->d(iz)==M->d()) {
      for (unsigned i=ie; i<M->e(iz)+ie; ++i) {
        switch(M->vz[iz].t) {
          case (FETRIANGLE):      dmesh.elmTypes[i] = ELM_SIMPLEX; dmesh.numElmFaces[i] = 3; dmesh.numElmNodes[i] = 3; break;
          case (FEQUADRILATERAL): dmesh.elmTypes[i] = ELM_QUAD;    dmesh.numElmFaces[i] = 4; dmesh.numElmNodes[i] = 4; break;
          case (FETETRAHEDRON):   dmesh.elmTypes[i] = ELM_SIMPLEX; dmesh.numElmFaces[i] = 4; dmesh.numElmNodes[i] = 4; break;
          case (FEBRICK):         dmesh.elmTypes[i] = ELM_HEX;     dmesh.numElmFaces[i] = 6; dmesh.numElmNodes[i] = 8; break;
          case (PRISM3):          dmesh.elmTypes[i] = ELM_PRISM;   dmesh.numElmFaces[i] = 5; dmesh.numElmNodes[i] = 5; break;
          case (PYRAMID4):        dmesh.elmTypes[i] = ELM_PYRAMID; dmesh.numElmFaces[i] = 5; dmesh.numElmNodes[i] = 5; break;
          default:                dmesh.elmTypes[i] = 0;
        }
        if (dmesh.elmTypes[i])
          dmesh.elmNorms[i].assign(dmesh.numElmFaces[i],std::vector< double >(M->d(),0.));
      }
      ie += M->e(iz);
    }
  }


  // set zones properties (if it's inner, if it's a wall)
  m_zprops.resize(M->z());
  for (unsigned i=0; i<M->z(); ++i)
    m_zprops[i].isinner = (M->d(i)==M->d());
  for (int i=0; i<x.nChildNode("wall"); ++i)
    m_zprops[ t_plas_aux::getzoneidx(x.getChildNode("wall",i).getAttribute< std::string >("zone"),*M) ].iswall = true;


  // mark boundary nodes
  {
#if 0
  std::vector< bool > nod_isbnd(M->n(),false);
  for (unsigned iz=0; iz<M->z(); ++iz) if (M->d(iz)==M->d()-1) {
    for (unsigned ie=0; ie<M->e(iz); ++ie) {
      std::vector< unsigned >::const_iterator in;
      for (in = M->vz[iz].e2n[ie].n.begin(); in!=M->vz[iz].e2n[ie].n.end(); ++in)
        nod_isbnd[*in] = true;
    }
  }


  m_zinner.resize(M->z());
  m_zbound.resize(M->z());
  for (unsigned i=0; i<M->z(); ++i) {
    if      (M->d(i)==M->d())   { m_zinner[i].nelem = M->e(i); m_zinner[i].etype = M->vz[i].t; }
    else if (M->d(i)==M->d()-1) { m_zbound[i].nelem = M->e(i); m_zbound[i].etype = M->vz[i].t; }
  }

  // mark boundary elements
  std::vector< std::vector< bool > > elm_isbnd(dmesh.numElm,std::vector< bool >(M->n(),false));
  for (unsigned iz=0, ib=0; iz<M->z(); ++iz) {
    if (M->d(iz)==M->d()-1) {
      for (unsigned jz=0, ie=0; jz<M->z(); ++jz)
        if (M->d(jz)==M->d()) {
          for (unsigned j=0; j<M->e(jz); ++j)
            for (std::vector< unsigned >::const_iterator n=M->vz[jz].e2n[j].n.begin(); n!=M->vz[jz].e2n[j].n.end(); ++n)
              elm_isbnd[ib][ie+j] = elm_isbnd[ib][ie+j] || nod_isbnd[ib][*n];
          ie += M->e(jz);
        }
      ++ib;
    }
  }
#endif
  }


  // set boundary elements
  dmesh.numBndFaces.assign(M->z(),0);
  dmesh.bndFaces   .resize(M->z());
  dmesh.bndDomElms .resize(M->z());
  for (unsigned iz=0, ib=0; iz<M->z(); ++iz) {
    if (M->d(iz)==M->d()-1) {

      // set number of boundary faces
      dmesh.numBndFaces[ib] = M->e(iz);

      // (iterate over this boundary's elements)
      typedef std::vector< unsigned > t_element_nodes;
      for (unsigned be=0; be<M->e(iz); ++be) {
// t_element_nodes &elm_bnd = M->vz[iz].e2n[be].n;
// FIXME ?
      }

      // set boundary faces "inner" element
      dmesh.bndDomElms[ib].assign(M->e(iz),0);
      for (unsigned jz=0, ie=0; jz<M->z(); ++jz) {
        if (M->d(jz)==M->d()) {
          ie += M->e(jz);
        }
      }

      // set boundary faces boundary face index (0-based)
      dmesh.bndFaces[ib].assign(M->e(iz),0);
      for (unsigned i=0; i<M->e(iz); ++i) {
// FIXME dmesh.bndFaces  [ib][i] = ?;
      }

      ++ib;
    }
  }
}


void t_plas::setFlowSolverParamOnInit(PLAS_FLOWSOLVER_PARAM *fp)
{
  fp->numDim       = M->d();
  fp->numUnk       = M->d()+2;
  fp->numNod       = M->n();
  fp->numElm       = dmesh.numElm;
  fp->numBnd       = M->z();
  fp->rhoCont      = dflow.rho;
  fp->muCont       = dflow.mu;
  fp->nuCont       = dflow.mu/dflow.rho;
  fp->cpCont       = dflow.cp;
  fp->kCont        = dflow.k;
  fp->dtEul        = dparam.dt;
  fp->minElmVolume = dmesh.minElmVolume;
  fp->maxElmVolume = dmesh.maxElmVolume;
  fp->writeOutput  = 1;

  fp->time = 0.;
  fp->iter = 0;
}


void t_plas::setFlowSolverParamOnTimeStep(PLAS_FLOWSOLVER_PARAM *fp)
{
  fp->time += dparam.dt;
  fp->iter =  dparam.iter;
}
