
#include <numeric>
#include "boost/progress.hpp"
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


  unsigned getvariableidx(const std::string& n, const m::mmesh& m)
  {
    for (unsigned r=0; r<m.v(); ++r)
      if (m.vn[r]==n)
        return r;
    std::cerr << "error: variable \"" << n << "\" not present!" << std::endl;
    throw 42;
    return 0;
  }


  unsigned getzoneidx(const std::string& n, const m::mmesh& m)
  {
    for (unsigned r=0; r<m.z(); ++r)
      if (m.vz[r].n==n)
        return r;
    std::cerr << "error: zone \"" << n << "\" not present!" << std::endl;
    throw 42;
    return 0;
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

  dparam.numIter  = dparam.numIter<0? 0      : dparam.numIter;
  dparam.dt       = dparam.dt<=0.?    1.e-12 : dparam.dt;
  plas::mdc = m::Create< PLAS_MATERIAL_DATA >(x.getAttribute< std::string >("continuum.material","water"));
  cout << "info: setup plas xml." << endl;


  cout << "info: setting quantities to provide..." << endl;

  // get variables to share
  m_quantity_idx.assign(ALL_QUANTITIES,-1);
  m_quantity_idx[COORD_X]    = t_plas_aux::getvariableidx("x",m);
  m_quantity_idx[COORD_Y]    = t_plas_aux::getvariableidx("y",m);
  m_quantity_idx[PRESSURE]   = t_plas_aux::getvariableidx("p",m);
  m_quantity_idx[VELOCITY_X] = t_plas_aux::getvariableidx("vx",m);
  m_quantity_idx[VELOCITY_Y] = t_plas_aux::getvariableidx("vy",m);
  if (m.d()>2) {
    m_quantity_idx[COORD_Z]    = t_plas_aux::getvariableidx("z",m);
    m_quantity_idx[VELOCITY_Z] = t_plas_aux::getvariableidx("vz",m);
  }

  // steady-state has time step n-1 quantities the same as time step n
  m_quantityold_idx = m_quantity_idx;

  cout << "info: setting quantities to provide." << endl;


  cout << "info: recreating data structures from m::mmesh..." << endl;
  M = &m;


  cout << "info: recreating: inner/boundary zones properties..." << endl;
  m_zinner_props.resize(M->z());
  m_zbound_props.resize(M->z());
  int
    ninnerelm = 0,  // total number of inner elements
    nboundelm = 0;  // total number of boundary elements
  for (unsigned iz=0; iz<M->z(); ++iz) {
    if (M->d(iz)==M->d()) {
      s_zoneprops &z = m_zinner_props[iz];
      z.nelems = M->e(iz);
      switch(M->vz[iz].t) {
      case (FETRIANGLE):      z.e_type = ELM_SIMPLEX; z.e_nfaces = 3; z.e_nnodes = 3; break;
      case (FEQUADRILATERAL): z.e_type = ELM_QUAD;    z.e_nfaces = 4; z.e_nnodes = 4; break;
      case (FETETRAHEDRON):   z.e_type = ELM_SIMPLEX; z.e_nfaces = 4; z.e_nnodes = 4; break;
      case (FEBRICK):         z.e_type = ELM_HEX;     z.e_nfaces = 6; z.e_nnodes = 8; break;
      case (PRISM3):          z.e_type = ELM_PRISM;   z.e_nfaces = 5; z.e_nnodes = 5; break;
      case (PYRAMID4):        z.e_type = ELM_PYRAMID; z.e_nfaces = 5; z.e_nnodes = 5; break;
      default: {}
      }
    }
    if (M->d(iz)==M->d()-1)
      m_zbound_props[iz].nelems = M->e(iz);

    ninnerelm += (M->d(iz)==M->d()?   M->e(iz):0);
    nboundelm += (M->d(iz)==M->d()-1? M->e(iz):0);
  }
  for (int i=0; i<x.nChildNode("wall"); ++i)
    m_zbound_props[ t_plas_aux::getzoneidx(x.getChildNode("wall",i).getAttribute< std::string >("zone"),*M) ].iswall = true;
  cout << "info: number of inner/boundary elements: " << ninnerelm << '/' << nboundelm << endl;
  cout << "info: recreating: inner/boundary zones properties." << endl;


  cout << "info: recreating: node-to-element connectivity..." << endl;
  dmesh.nodElms.assign(M->n(),std::vector< int >());
  for (size_t iz=0, ielm=0; iz<m_zinner_props.size(); ++iz)
    for (int ie=0; ie<m_zinner_props[iz].nelems; ++ie, ++ielm)
      for (std::vector< unsigned >::const_iterator n=M->vz[iz].e2n[ie].n.begin(); n!=M->vz[iz].e2n[ie].n.end(); ++n)
        dmesh.nodElms[ *n ].push_back(ielm);
  cout << "info: recreating: node-to-element connectivity." << endl;


  cout << "info: recreating: element-to-element (sharing a face) connectivity..." << endl;
  dmesh.elmNeighbs.resize(ninnerelm);
  std::vector< int >
    ifacenodes(4,-1),
    jfacenodes(4,-1);
  std::vector< std::vector< int > > v_isboudaryelm(M->z());
  boost::progress_display pbar(ninnerelm);
  for (size_t iz=0, ielm=0; iz<m_zinner_props.size(); ++iz) {
    if (m_zinner_props[iz].nelems)
      v_isboudaryelm[iz].assign(m_zinner_props[iz].nelems,0);
    for (int ie=0; ie<m_zinner_props[iz].nelems; ++ie, ++ielm, ++pbar) {
      dmesh.elmNeighbs[ielm].assign(m_zinner_props[iz].e_nfaces,-1);
      for (int ifac=0; ifac<m_zinner_props[iz].e_nfaces; ++ifac) {

        // get nodes in the face
        plasdriver_GetFaceNodes(iz,ie,ifac,&ifacenodes[0]);
        std::sort(ifacenodes.begin(),ifacenodes.end());

        // get list of searcheable elements (sharing a node), absolute index
        std::vector< int > v_ielm;
        for (std::vector< int >::const_iterator n=ifacenodes.begin(); n!=ifacenodes.end(); ++n)
          if (*n!=-1)
            v_ielm.insert(v_ielm.end(),dmesh.nodElms[*n].begin(),dmesh.nodElms[*n].end());
        sort(v_ielm.begin(),v_ielm.end());
        v_ielm.erase(unique(v_ielm.begin(),v_ielm.end()),v_ielm.end());

        // convert absolute indices to relative (zone and internal elem. index)
        std::vector< int >
          v_jelem_z,
          v_jelem_e,
          v_jelem_elm;
        for (std::vector< int >::const_iterator elm=v_ielm.begin(); elm!=v_ielm.end(); ++elm)
          for (int jz=0, nelems=0; jz<(int) m_zinner_props.size(); ++jz) {
            nelems += m_zinner_props[jz].nelems;
            if (nelems > *elm) {
              v_jelem_z.push_back(jz);
              v_jelem_e.push_back(*elm - nelems + m_zinner_props[jz].nelems);
              v_jelem_elm.push_back(*elm);
              break;
            }
          }
        v_ielm.clear();

        // find the element with a face with the same nodes but different index
        for (size_t j=0; j<v_jelem_z.size(); ++j) {
          const int
            _z = v_jelem_z[j],
            _e = v_jelem_e[j];
          size_t jelm = v_jelem_elm[j];
          for (int jfac=0; jfac<m_zinner_props[_z].e_nfaces && dmesh.elmNeighbs[ielm][ifac]<0; ++jfac) {
            plasdriver_GetFaceNodes(_z,_e,jfac,&jfacenodes[0]);
            std::sort(jfacenodes.begin(),jfacenodes.end());
            dmesh.elmNeighbs[ielm][ifac] = (ielm!=jelm && ifacenodes==jfacenodes? jelm : -1);
          }
        }

        v_isboudaryelm[iz][ie] += (dmesh.elmNeighbs[ielm][ifac]<0? 1:0);
      }
    }
  }
  cout << "info: recreating: element-to-element (sharing a face) connectivity." << endl;


  cout << "info: recreating: boundary elements (to inner elements)..." << endl;
  dmesh.bndFaces  .resize(M->z());
  dmesh.bndDomElms.resize(M->z());
  pbar.restart(nboundelm);
  for (unsigned iz=0, ielm=0; iz<m_zbound_props.size(); ++iz) {
    if (m_zbound_props[iz].nelems) {
      dmesh.bndFaces  [iz].assign(m_zbound_props[iz].nelems,-1);
      dmesh.bndDomElms[iz].assign(m_zbound_props[iz].nelems,-1);
    }
    for (int ie=0; ie<m_zbound_props[iz].nelems; ++ie, ++ielm, ++pbar) {

      // create digestible face nodes
      std::vector< int > ifacenodes(4,-1);
      for (unsigned i=0; i<M->vz[iz].e2n[ie].n.size() && i<4; ++i)
        ifacenodes[i] = M->vz[iz].e2n[ie].n[i];
      std::sort(ifacenodes.begin(),ifacenodes.end());

      // create faces to search for
      std::vector< int > jfacenodes(4);
      for (size_t jz=0, jelm=0; jz<m_zinner_props.size(); ++jz) {
        for (int je=0; je<m_zinner_props[jz].nelems; ++je, ++jelm) {
          if (v_isboudaryelm[jz][je]) {
            for (int jfac=0; jfac<m_zinner_props[jz].e_nfaces; ++jfac) {
              plasdriver_GetFaceNodes(jz,je,jfac,&jfacenodes[0]);
              std::sort(jfacenodes.begin(),jfacenodes.end());

              if (ifacenodes==jfacenodes) {
                dmesh.bndFaces  [iz][ie] = jfac;
                dmesh.bndDomElms[iz][ie] = jelm;
                --v_isboudaryelm[jz][je];
                break;
              }

            }
          }
        }
      }

    }
  }
  cout << "info: recreating: boundary inner elements (to inner elements)." << endl;


  cout << "info: recreating: inner element normals..." << endl;

  // allocate by number of (inner) elements, (inner) elements faces
  dmesh.elmNorms.resize(ninnerelm);
  for (size_t iz=0, e=0; iz<m_zinner_props.size(); ++iz)
    for (int ie=0; ie<m_zinner_props[iz].nelems; ++ie, ++e)
      dmesh.elmNorms[e].assign( m_zinner_props[iz].e_nfaces, std::vector< double >(M->d(),0.) );

  // calculate face normals
  int fnodes[4];
  for (size_t iz=0, e=0; iz<m_zinner_props.size(); ++iz) {
    for (int ie=0; ie<m_zinner_props[iz].nelems; ++ie, ++e) {
      for (int f=0; f<m_zinner_props[iz].e_nfaces; ++f) {

        fnodes[0] = fnodes[1] = fnodes[2] = fnodes[3] = -1;
        plasdriver_GetFaceNodes(iz,ie,f,fnodes);
        const int e_type = m_zinner_props[iz].e_type;

        if (e_type==ELM_SIMPLEX && M->d()==2) {
          dmesh.elmNorms[e][f][0] = M->vv[1][fnodes[0]] - M->vv[1][fnodes[1]];
          dmesh.elmNorms[e][f][1] = M->vv[0][fnodes[1]] - M->vv[0][fnodes[0]];
        }
        else if ((e_type==ELM_SIMPLEX && M->d()==3)
              || (e_type==ELM_PRISM   && f>2)
              || (e_type==ELM_PYRAMID && f>0)) {
          dmesh.elmNorms[e][f][0] =
            0.5*((M->vv[1][fnodes[2]] - M->vv[1][fnodes[0]])
                *(M->vv[2][fnodes[1]] - M->vv[2][fnodes[0]])
                -(M->vv[1][fnodes[1]] - M->vv[1][fnodes[0]])
                *(M->vv[2][fnodes[2]] - M->vv[2][fnodes[0]]));
          dmesh.elmNorms[e][f][1] =
            0.5*((M->vv[2][fnodes[2]] - M->vv[2][fnodes[0]])
                *(M->vv[0][fnodes[1]] - M->vv[0][fnodes[0]])
                -(M->vv[2][fnodes[1]] - M->vv[2][fnodes[0]])
                *(M->vv[0][fnodes[2]] - M->vv[0][fnodes[0]]));
          dmesh.elmNorms[e][f][2] =
            0.5*((M->vv[0][fnodes[2]] - M->vv[0][fnodes[0]])
                *(M->vv[1][fnodes[1]] - M->vv[1][fnodes[0]])
                -(M->vv[0][fnodes[1]] - M->vv[0][fnodes[0]])
                *(M->vv[1][fnodes[2]] - M->vv[1][fnodes[0]]));
        }
        else if(e_type==ELM_QUAD) {
          dmesh.elmNorms[e][f][0] = 2.0*(M->vv[1][fnodes[0]] - M->vv[1][fnodes[1]]);
          dmesh.elmNorms[e][f][1] = 2.0*(M->vv[0][fnodes[1]] - M->vv[0][fnodes[0]]);
        }
        else if (e_type==ELM_HEX
             || (e_type==ELM_PRISM   && f<=2)
             || (e_type==ELM_PYRAMID && f==0)) {
          dmesh.elmNorms[e][f][0] =
              ((M->vv[1][fnodes[3]] - M->vv[1][fnodes[0]])
              *(M->vv[2][fnodes[1]] - M->vv[2][fnodes[0]])
              -(M->vv[1][fnodes[1]] - M->vv[1][fnodes[0]])
              *(M->vv[2][fnodes[3]] - M->vv[2][fnodes[0]]));
          dmesh.elmNorms[e][f][1] =
              ((M->vv[2][fnodes[3]] - M->vv[2][fnodes[0]])
              *(M->vv[0][fnodes[1]] - M->vv[0][fnodes[0]])
              -(M->vv[2][fnodes[1]] - M->vv[2][fnodes[0]])
              *(M->vv[0][fnodes[3]] - M->vv[0][fnodes[0]]));
          dmesh.elmNorms[e][f][2] =
              ((M->vv[0][fnodes[3]] - M->vv[0][fnodes[0]])
              *(M->vv[1][fnodes[1]] - M->vv[1][fnodes[0]])
              -(M->vv[0][fnodes[1]] - M->vv[0][fnodes[0]])
              *(M->vv[1][fnodes[3]] - M->vv[1][fnodes[0]]));
        }

      }
    }
  }
  cout << "info: recreating: inner element normals." << endl;


  cout << "info: recreating: inner element volumes..." << endl;
  dmesh.elmVolumes.assign(ninnerelm,0.);
  double c2[3][2], c3[4][3];
  unsigned ielm = 0;
  for (size_t iz=0; iz<m_zinner_props.size(); ++iz) {
    for (int ie=0; ie<m_zinner_props[iz].nelems; ++ie) {
      double &volume = dmesh.elmVolumes[ielm+ie];

      if (m_zinner_props[iz].e_type==ELM_SIMPLEX && M->d()==2){

        c2[0][0] = M->vv[0][ M->vz[iz].e2n[ie].n[0] ];
        c2[0][1] = M->vv[1][ M->vz[iz].e2n[ie].n[0] ];
        c2[1][0] = M->vv[0][ M->vz[iz].e2n[ie].n[1] ];
        c2[1][1] = M->vv[1][ M->vz[iz].e2n[ie].n[1] ];
        c2[2][0] = M->vv[0][ M->vz[iz].e2n[ie].n[2] ];
        c2[2][1] = M->vv[1][ M->vz[iz].e2n[ie].n[2] ];
        volume = plasdriver_CalcAreaTriangle(c2);

      }
      else if (m_zinner_props[iz].e_type==ELM_SIMPLEX && M->d()==3){

        c3[0][0] = M->vv[0][ M->vz[iz].e2n[ie].n[0] ];
        c3[0][1] = M->vv[1][ M->vz[iz].e2n[ie].n[0] ];
        c3[0][2] = M->vv[2][ M->vz[iz].e2n[ie].n[0] ];
        c3[1][0] = M->vv[0][ M->vz[iz].e2n[ie].n[1] ];
        c3[1][1] = M->vv[1][ M->vz[iz].e2n[ie].n[1] ];
        c3[1][2] = M->vv[2][ M->vz[iz].e2n[ie].n[1] ];
        c3[2][0] = M->vv[0][ M->vz[iz].e2n[ie].n[2] ];
        c3[2][1] = M->vv[1][ M->vz[iz].e2n[ie].n[2] ];
        c3[2][2] = M->vv[2][ M->vz[iz].e2n[ie].n[2] ];
        c3[3][0] = M->vv[0][ M->vz[iz].e2n[ie].n[3] ];
        c3[3][1] = M->vv[1][ M->vz[iz].e2n[ie].n[3] ];
        c3[3][2] = M->vv[2][ M->vz[iz].e2n[ie].n[3] ];
        volume = plasdriver_CalcVolumeTetra(c3);

      }
      else if (m_zinner_props[iz].e_type==ELM_QUAD){

        c2[0][0] = M->vv[0][ M->vz[iz].e2n[ie].n[0] ];
        c2[0][1] = M->vv[1][ M->vz[iz].e2n[ie].n[0] ];
        c2[1][0] = M->vv[0][ M->vz[iz].e2n[ie].n[1] ];
        c2[1][1] = M->vv[1][ M->vz[iz].e2n[ie].n[1] ];
        c2[2][0] = M->vv[0][ M->vz[iz].e2n[ie].n[2] ];
        c2[2][1] = M->vv[1][ M->vz[iz].e2n[ie].n[2] ];
        volume = plasdriver_CalcAreaTriangle(c2);

        c2[0][0] = M->vv[0][ M->vz[iz].e2n[ie].n[0] ];
        c2[0][1] = M->vv[1][ M->vz[iz].e2n[ie].n[0] ];
        c2[1][0] = M->vv[0][ M->vz[iz].e2n[ie].n[2] ];
        c2[1][1] = M->vv[1][ M->vz[iz].e2n[ie].n[2] ];
        c2[2][0] = M->vv[0][ M->vz[iz].e2n[ie].n[3] ];
        c2[2][1] = M->vv[1][ M->vz[iz].e2n[ie].n[3] ];
        volume += plasdriver_CalcAreaTriangle(c2);

      }
      else if (m_zinner_props[iz].e_type==ELM_HEX){

        c3[0][0] = M->vv[0][ M->vz[iz].e2n[ie].n[0] ];
        c3[0][1] = M->vv[1][ M->vz[iz].e2n[ie].n[0] ];
        c3[0][2] = M->vv[2][ M->vz[iz].e2n[ie].n[0] ];
        c3[1][0] = M->vv[0][ M->vz[iz].e2n[ie].n[1] ];
        c3[1][1] = M->vv[1][ M->vz[iz].e2n[ie].n[1] ];
        c3[1][2] = M->vv[2][ M->vz[iz].e2n[ie].n[1] ];
        c3[2][0] = M->vv[0][ M->vz[iz].e2n[ie].n[3] ];
        c3[2][1] = M->vv[1][ M->vz[iz].e2n[ie].n[3] ];
        c3[2][2] = M->vv[2][ M->vz[iz].e2n[ie].n[3] ];
        c3[3][0] = M->vv[0][ M->vz[iz].e2n[ie].n[5] ];
        c3[3][1] = M->vv[1][ M->vz[iz].e2n[ie].n[5] ];
        c3[3][2] = M->vv[2][ M->vz[iz].e2n[ie].n[5] ];
        volume = plasdriver_CalcVolumeTetra(c3);

        c3[0][0] = M->vv[0][ M->vz[iz].e2n[ie].n[0] ];
        c3[0][1] = M->vv[1][ M->vz[iz].e2n[ie].n[0] ];
        c3[0][2] = M->vv[2][ M->vz[iz].e2n[ie].n[0] ];
        c3[1][0] = M->vv[0][ M->vz[iz].e2n[ie].n[3] ];
        c3[1][1] = M->vv[1][ M->vz[iz].e2n[ie].n[3] ];
        c3[1][2] = M->vv[2][ M->vz[iz].e2n[ie].n[3] ];
        c3[2][0] = M->vv[0][ M->vz[iz].e2n[ie].n[6] ];
        c3[2][1] = M->vv[1][ M->vz[iz].e2n[ie].n[6] ];
        c3[2][2] = M->vv[2][ M->vz[iz].e2n[ie].n[6] ];
        c3[3][0] = M->vv[0][ M->vz[iz].e2n[ie].n[5] ];
        c3[3][1] = M->vv[1][ M->vz[iz].e2n[ie].n[5] ];
        c3[3][2] = M->vv[2][ M->vz[iz].e2n[ie].n[5] ];
        volume += plasdriver_CalcVolumeTetra(c3);

        c3[0][0] = M->vv[0][ M->vz[iz].e2n[ie].n[0] ];
        c3[0][1] = M->vv[1][ M->vz[iz].e2n[ie].n[0] ];
        c3[0][2] = M->vv[2][ M->vz[iz].e2n[ie].n[0] ];
        c3[1][0] = M->vv[0][ M->vz[iz].e2n[ie].n[3] ];
        c3[1][1] = M->vv[1][ M->vz[iz].e2n[ie].n[3] ];
        c3[1][2] = M->vv[2][ M->vz[iz].e2n[ie].n[3] ];
        c3[2][0] = M->vv[0][ M->vz[iz].e2n[ie].n[2] ];
        c3[2][1] = M->vv[1][ M->vz[iz].e2n[ie].n[2] ];
        c3[2][2] = M->vv[2][ M->vz[iz].e2n[ie].n[2] ];
        c3[3][0] = M->vv[0][ M->vz[iz].e2n[ie].n[6] ];
        c3[3][1] = M->vv[1][ M->vz[iz].e2n[ie].n[6] ];
        c3[3][2] = M->vv[2][ M->vz[iz].e2n[ie].n[6] ];
        volume += plasdriver_CalcVolumeTetra(c3);

        c3[0][0] = M->vv[0][ M->vz[iz].e2n[ie].n[0] ];
        c3[0][1] = M->vv[1][ M->vz[iz].e2n[ie].n[0] ];
        c3[0][2] = M->vv[2][ M->vz[iz].e2n[ie].n[0] ];
        c3[1][0] = M->vv[0][ M->vz[iz].e2n[ie].n[5] ];
        c3[1][1] = M->vv[1][ M->vz[iz].e2n[ie].n[5] ];
        c3[1][2] = M->vv[2][ M->vz[iz].e2n[ie].n[5] ];
        c3[2][0] = M->vv[0][ M->vz[iz].e2n[ie].n[6] ];
        c3[2][1] = M->vv[1][ M->vz[iz].e2n[ie].n[6] ];
        c3[2][2] = M->vv[2][ M->vz[iz].e2n[ie].n[6] ];
        c3[3][0] = M->vv[0][ M->vz[iz].e2n[ie].n[4] ];
        c3[3][1] = M->vv[1][ M->vz[iz].e2n[ie].n[4] ];
        c3[3][2] = M->vv[2][ M->vz[iz].e2n[ie].n[4] ];
        volume += plasdriver_CalcVolumeTetra(c3);

        c3[0][0] = M->vv[0][ M->vz[iz].e2n[ie].n[3] ];
        c3[0][1] = M->vv[1][ M->vz[iz].e2n[ie].n[3] ];
        c3[0][2] = M->vv[2][ M->vz[iz].e2n[ie].n[3] ];
        c3[1][0] = M->vv[0][ M->vz[iz].e2n[ie].n[6] ];
        c3[1][1] = M->vv[1][ M->vz[iz].e2n[ie].n[6] ];
        c3[1][2] = M->vv[2][ M->vz[iz].e2n[ie].n[6] ];
        c3[2][0] = M->vv[0][ M->vz[iz].e2n[ie].n[5] ];
        c3[2][1] = M->vv[1][ M->vz[iz].e2n[ie].n[5] ];
        c3[2][2] = M->vv[2][ M->vz[iz].e2n[ie].n[5] ];
        c3[3][0] = M->vv[0][ M->vz[iz].e2n[ie].n[7] ];
        c3[3][1] = M->vv[1][ M->vz[iz].e2n[ie].n[7] ];
        c3[3][2] = M->vv[2][ M->vz[iz].e2n[ie].n[7] ];
        volume += plasdriver_CalcVolumeTetra(c3);

      }
      else if (m_zinner_props[iz].e_type==ELM_PRISM){

        c3[0][0] = M->vv[0][ M->vz[iz].e2n[ie].n[0] ];
        c3[0][1] = M->vv[1][ M->vz[iz].e2n[ie].n[0] ];
        c3[0][2] = M->vv[2][ M->vz[iz].e2n[ie].n[0] ];
        c3[1][0] = M->vv[0][ M->vz[iz].e2n[ie].n[1] ];
        c3[1][1] = M->vv[1][ M->vz[iz].e2n[ie].n[1] ];
        c3[1][2] = M->vv[2][ M->vz[iz].e2n[ie].n[1] ];
        c3[2][0] = M->vv[0][ M->vz[iz].e2n[ie].n[2] ];
        c3[2][1] = M->vv[1][ M->vz[iz].e2n[ie].n[2] ];
        c3[2][2] = M->vv[2][ M->vz[iz].e2n[ie].n[2] ];
        c3[3][0] = M->vv[0][ M->vz[iz].e2n[ie].n[5] ];
        c3[3][1] = M->vv[1][ M->vz[iz].e2n[ie].n[5] ];
        c3[3][2] = M->vv[2][ M->vz[iz].e2n[ie].n[5] ];
        volume = plasdriver_CalcVolumeTetra(c3);

        c3[0][0] = M->vv[0][ M->vz[iz].e2n[ie].n[0] ];
        c3[0][1] = M->vv[1][ M->vz[iz].e2n[ie].n[0] ];
        c3[0][2] = M->vv[2][ M->vz[iz].e2n[ie].n[0] ];
        c3[1][0] = M->vv[0][ M->vz[iz].e2n[ie].n[1] ];
        c3[1][1] = M->vv[1][ M->vz[iz].e2n[ie].n[1] ];
        c3[1][2] = M->vv[2][ M->vz[iz].e2n[ie].n[1] ];
        c3[2][0] = M->vv[0][ M->vz[iz].e2n[ie].n[5] ];
        c3[2][1] = M->vv[1][ M->vz[iz].e2n[ie].n[5] ];
        c3[2][2] = M->vv[2][ M->vz[iz].e2n[ie].n[5] ];
        c3[3][0] = M->vv[0][ M->vz[iz].e2n[ie].n[4] ];
        c3[3][1] = M->vv[1][ M->vz[iz].e2n[ie].n[4] ];
        c3[3][2] = M->vv[2][ M->vz[iz].e2n[ie].n[4] ];
        volume += plasdriver_CalcVolumeTetra(c3);

        c3[0][0] = M->vv[0][ M->vz[iz].e2n[ie].n[0] ];
        c3[0][1] = M->vv[1][ M->vz[iz].e2n[ie].n[0] ];
        c3[0][2] = M->vv[2][ M->vz[iz].e2n[ie].n[0] ];
        c3[1][0] = M->vv[0][ M->vz[iz].e2n[ie].n[4] ];
        c3[1][1] = M->vv[1][ M->vz[iz].e2n[ie].n[4] ];
        c3[1][2] = M->vv[2][ M->vz[iz].e2n[ie].n[4] ];
        c3[2][0] = M->vv[0][ M->vz[iz].e2n[ie].n[5] ];
        c3[2][1] = M->vv[1][ M->vz[iz].e2n[ie].n[5] ];
        c3[2][2] = M->vv[2][ M->vz[iz].e2n[ie].n[5] ];
        c3[3][0] = M->vv[0][ M->vz[iz].e2n[ie].n[3] ];
        c3[3][1] = M->vv[1][ M->vz[iz].e2n[ie].n[3] ];
        c3[3][2] = M->vv[2][ M->vz[iz].e2n[ie].n[3] ];
        volume += plasdriver_CalcVolumeTetra(c3);

      }
      else if (m_zinner_props[iz].e_type==ELM_PYRAMID){

        c3[0][0] = M->vv[0][ M->vz[iz].e2n[ie].n[0] ];
        c3[0][1] = M->vv[1][ M->vz[iz].e2n[ie].n[0] ];
        c3[0][2] = M->vv[2][ M->vz[iz].e2n[ie].n[0] ];
        c3[1][0] = M->vv[0][ M->vz[iz].e2n[ie].n[1] ];
        c3[1][1] = M->vv[1][ M->vz[iz].e2n[ie].n[1] ];
        c3[1][2] = M->vv[2][ M->vz[iz].e2n[ie].n[1] ];
        c3[2][0] = M->vv[0][ M->vz[iz].e2n[ie].n[3] ];
        c3[2][1] = M->vv[1][ M->vz[iz].e2n[ie].n[3] ];
        c3[2][2] = M->vv[2][ M->vz[iz].e2n[ie].n[3] ];
        c3[3][0] = M->vv[0][ M->vz[iz].e2n[ie].n[4] ];
        c3[3][1] = M->vv[1][ M->vz[iz].e2n[ie].n[4] ];
        c3[3][2] = M->vv[2][ M->vz[iz].e2n[ie].n[4] ];
        volume = plasdriver_CalcVolumeTetra(c3);

        c3[0][0] = M->vv[0][ M->vz[iz].e2n[ie].n[0] ];
        c3[0][1] = M->vv[1][ M->vz[iz].e2n[ie].n[0] ];
        c3[0][2] = M->vv[2][ M->vz[iz].e2n[ie].n[0] ];
        c3[1][0] = M->vv[0][ M->vz[iz].e2n[ie].n[3] ];
        c3[1][1] = M->vv[1][ M->vz[iz].e2n[ie].n[3] ];
        c3[1][2] = M->vv[2][ M->vz[iz].e2n[ie].n[3] ];
        c3[2][0] = M->vv[0][ M->vz[iz].e2n[ie].n[2] ];
        c3[2][1] = M->vv[1][ M->vz[iz].e2n[ie].n[2] ];
        c3[2][2] = M->vv[2][ M->vz[iz].e2n[ie].n[2] ];
        c3[3][0] = M->vv[0][ M->vz[iz].e2n[ie].n[4] ];
        c3[3][1] = M->vv[1][ M->vz[iz].e2n[ie].n[4] ];
        c3[3][2] = M->vv[2][ M->vz[iz].e2n[ie].n[4] ];
        volume += plasdriver_CalcVolumeTetra(c3);

      }

    }
    ielm += m_zinner_props[iz].nelems;
  }
  cout << "info: recreating: inner element volumes." << endl;


  cout << "info: recreating: nodal volumes (dual cell)..." << endl;
  dmesh.nodVolumes.assign(M->n(),0.);
  for (unsigned n=0; n<M->n(); ++n)
    for (size_t e=0; e<dmesh.nodElms[n].size(); ++e)
      dmesh.nodVolumes[n] += dmesh.elmVolumes[dmesh.nodElms[n][e]]/dmesh.nodElms[n].size();
  cout << "info: recreating: nodal volumes (dual cell)." << endl;


  cout << "info: recreating data structures from m::mmesh." << endl;


  screenOutput("initializing PLaS...");
  plas::initialize(x);
  screenOutput("initializing PLaS.");


  screenOutput("perform PLaS iterations...");
  for (dparam.iter=1; dparam.iter<=dparam.numIter; ++dparam.iter)
    plas::run();
  screenOutput("perform PLaS iterations.");
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
void t_plas::plasdriver_GetFaceNodes(int iz, int ie, int face, int *nodes)
{
  nodes[0] = nodes[1] = nodes[2] = nodes[3] = -1;
  const std::vector< unsigned > &en = M->vz[iz].e2n[ie].n;

  switch (m_zinner_props[iz].e_type) {
  case ELM_SIMPLEX:
    if (M->d()==2) {
      switch (face) {
      case 0: { nodes[0]=en[0]; nodes[1]=en[1]; break; }
      case 1: { nodes[0]=en[1]; nodes[1]=en[2]; break; }
      case 2: { nodes[0]=en[2]; nodes[1]=en[0]; break; }
      }
    }
    else if (M->d()==3) {
      switch (face) {
      case 0: { nodes[0]=en[1]; nodes[1]=en[0]; nodes[2]=en[2]; break; }
      case 1: { nodes[0]=en[0]; nodes[1]=en[1]; nodes[2]=en[3]; break; }
      case 2: { nodes[0]=en[1]; nodes[1]=en[2]; nodes[2]=en[3]; break; }
      case 3: { nodes[0]=en[2]; nodes[1]=en[0]; nodes[2]=en[3]; break; }
      }
    } break;
  case ELM_QUAD:
    switch (face) {
    case 0: { nodes[0]=en[0]; nodes[1]=en[1]; break; }
    case 1: { nodes[0]=en[1]; nodes[1]=en[2]; break; }
    case 2: { nodes[0]=en[2]; nodes[1]=en[3]; break; }
    case 3: { nodes[0]=en[3]; nodes[1]=en[0]; break; }
    } break;
  case ELM_HEX:
    switch (face) {
    case 0: { nodes[0]=en[0]; nodes[1]=en[1]; nodes[2]=en[5]; nodes[3]=en[4]; break; }
    case 1: { nodes[0]=en[1]; nodes[1]=en[3]; nodes[2]=en[7]; nodes[3]=en[5]; break; }
    case 2: { nodes[0]=en[3]; nodes[1]=en[2]; nodes[2]=en[6]; nodes[3]=en[7]; break; }
    case 3: { nodes[0]=en[2]; nodes[1]=en[0]; nodes[2]=en[4]; nodes[3]=en[6]; break; }
    case 4: { nodes[0]=en[1]; nodes[1]=en[0]; nodes[2]=en[2]; nodes[3]=en[3]; break; }
    case 5: { nodes[0]=en[4]; nodes[1]=en[5]; nodes[2]=en[7]; nodes[3]=en[6]; break; }
    } break;
  case ELM_PRISM:
    switch (face) {
    case 0: { nodes[0]=en[0]; nodes[1]=en[1]; nodes[2]=en[4]; nodes[3]=en[3]; break; }
    case 1: { nodes[0]=en[1]; nodes[1]=en[2]; nodes[2]=en[5]; nodes[3]=en[4]; break; }
    case 2: { nodes[0]=en[2]; nodes[1]=en[0]; nodes[2]=en[3]; nodes[3]=en[5]; break; }
    case 3: { nodes[0]=en[0]; nodes[1]=en[2]; nodes[2]=en[1]; break; }
    case 4: { nodes[0]=en[3]; nodes[1]=en[4]; nodes[2]=en[5]; break; }
    } break;
  case ELM_PYRAMID:
    switch (face) {
    case 0: { nodes[0]=en[0]; nodes[1]=en[2]; nodes[2]=en[3]; nodes[3]=en[1]; break; }
    case 1: { nodes[0]=en[0]; nodes[1]=en[1]; nodes[2]=en[4]; break; }
    case 2: { nodes[0]=en[1]; nodes[1]=en[3]; nodes[2]=en[4]; break; }
    case 3: { nodes[0]=en[3]; nodes[1]=en[2]; nodes[2]=en[4]; break; }
    case 4: { nodes[0]=en[2]; nodes[1]=en[0]; nodes[2]=en[4]; break; }
    } break;
  }
}


void t_plas::setFlowSolverParamOnInit(PLAS_FLOWSOLVER_PARAM *fp)
{
  int numElm = 0;
  for (size_t i=0; i<m_zinner_props.size(); ++i)
    numElm += m_zinner_props[i].nelems;

  fp->numDim = M->d();
  fp->numUnk = M->d()+2;
  fp->numNod = M->n();
  fp->numElm = numElm;
  fp->numBnd = M->z();

  fp->dtEul        = dparam.dt;
  fp->minElmVolume = *std::min_element(dmesh.elmVolumes.begin(),dmesh.elmVolumes.end());
  fp->maxElmVolume = *std::max_element(dmesh.elmVolumes.begin(),dmesh.elmVolumes.end());

  fp->time = 0.;
  fp->iter = 0;
}


void t_plas::setFlowSolverParamOnTimeStep(PLAS_FLOWSOLVER_PARAM *fp)
{
  fp->time += dparam.dt;
  fp->iter  = dparam.iter;
}
