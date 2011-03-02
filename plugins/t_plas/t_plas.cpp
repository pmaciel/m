
#include <numeric>
#include "boost/progress.hpp"
#include "ext/xmlParser.h"
#include "mfactory.h"
#include "t_plas.h"


m::Register< m::mtransform,t_plas > mt_plas( 6,
  "-tplas", "[str] filename or string with xml formatted as:",
  "", "<plas",
  "", " iterations=\"[int]\" (default 1)",
  "", " dt=\"[real]\" (default 1.)",
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


void t_plas::transform(GetPot& o, m::mmesh& m)
{
  screenOutput("setup plas xml...");
  const std::string o_xml = o.get(o.inc_cursor(),"");
  XMLNode x = ((o_xml.size()? o_xml[0]:'?')=='<')? XMLNode::parseString(o_xml.c_str(),"plas")
                                                 : XMLNode::openFileHelper(o_xml.c_str(),"plas");
  m_numiter = std::max(0,     x.getAttribute< int    >("iterations",1 ));
  m_dt      = std::max(1.e-12,x.getAttribute< double >("dt",        1.));
  m_ziswall.assign(m.z(),false);
  for (int i=0; i<x.nChildNode("wall"); ++i)
    m_ziswall[ t_plas_aux::getzoneidx(x.getChildNode("wall",i).getAttribute< std::string >("zone"),m) ] = true;
  screenOutput("setup plas xml.");


  screenOutput("setting quantities to provide...");
  // get variables to share
  m_quantity_idx.assign(ALL_QUANTITIES,-1);
  m_quantity_idx[COORD_X]      = t_plas_aux::getvariableidx("x",m);
  m_quantity_idx[COORD_Y]      = t_plas_aux::getvariableidx("y",m);
  m_quantity_idx[PRESSURE]     = t_plas_aux::getvariableidx("p",m);
  m_quantity_idx[VELOCITY_X]   = t_plas_aux::getvariableidx("vx",m);
  m_quantity_idx[VELOCITY_Y]   = t_plas_aux::getvariableidx("vy",m);
  if (m.d()>2) {
    m_quantity_idx[COORD_Z]    = t_plas_aux::getvariableidx("z",m);
    m_quantity_idx[VELOCITY_Z] = t_plas_aux::getvariableidx("vz",m);
  }
  screenOutput("setting quantities to provide.");


  // share the data structure accross member functions
  M = &m;


  screenOutput("initializing PLaS...");
  plas::initialize(x);
  screenOutput("initializing PLaS.");


  screenOutput("perform PLaS iterations...");
  for (m_iter=1; m_iter<=m_numiter; ++m_iter)
    plas::run();
  screenOutput("perform PLaS iterations.");
}


void t_plas::setFlowSolverParamOnInit(PLAS_FLOWSOLVER_PARAM *fp)
{
  fp->nInnerElements.assign(M->z(),0);
  fp->nBoundElements.assign(M->z(),0);
  for (size_t i=0; i<M->z(); ++i) {
    fp->nInnerElements[i] = (M->d()==M->d(i)?   M->e(i) : 0);
    fp->nBoundElements[i] = (M->d()==M->d(i)-1? M->e(i) : 0);
  }
  fp->numDim = M->d();
  fp->numUnk = M->d()+2;
  fp->numNod = M->n();
  fp->dtEul  = m_dt;

  fp->time = 0.;
  fp->iter = 0;
}
