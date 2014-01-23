
#include "ext/xmlParser.h"
#include "mfactory.h"
#include "t_pillaz.h"


m::Register< m::mtransform,t_pillaz > mt_pillaz( 6,
  "-tpillaz", "[str] filename or string with xml formatted as:",
  "", "<pillaz",
  "", " iterations=\"[int]\" (default 1)",
  "", " dt=\"[real]\" (default 1.)",
  "", " <wall zone=\"[str]\"/> (walls, 1 or more)",
  "", "</pillaz>" );


namespace t_pillaz_aux {


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


void t_pillaz::transform(GetPot& o, m::mmesh& m, const XMLNode& x, const XMLNode& x)
{

  // setup xml options
  const std::string o_xml = o.get(o.inc_cursor(),"");
  XMLNode x = ((o_xml.size()? o_xml[0]:'?')=='<')? XMLNode::parseString(o_xml.c_str(),"pillaz")
                                                 : XMLNode::openFileHelper(o_xml.c_str(),"pillaz");

  // set number of iterations and time-step
  m_numiter = std::max(0,     x.getAttribute< int    >("iterations",1 ));
  m_dt      = std::max(1.e-12,x.getAttribute< double >("dt",        1.));

  // setup walls
  m_ziswall.assign(m.z(),false);
  for (int i=0; i<x.nChildNode("wall"); ++i)
    m_ziswall[ t_pillaz_aux::getzoneidx(x.getChildNode("wall",i).getAttribute< std::string >("zone"),m) ] = true;

  // setup data structure (accross member functions)
  M = &m;

  // setup quantities to provide
  m_quantityscalar.assign(ALL_QSCALAR,-1);
  m_quantityvector.assign(ALL_QVECTOR,-1);
  m_quantityscalar[PRESSURE] = t_pillaz_aux::getvariableidx("p", m);
  m_quantityvector[VELOCITY] = t_pillaz_aux::getvariableidx("vx",m);
  m_quantityvector[COORD]    = t_pillaz_aux::getvariableidx("x", m);

  // run
  pillaz::initialize(x);
  for (int i=1; i<=m_numiter; ++i)
    pillaz::run();
}


void t_pillaz::setFlowSolverParamOnInit(PILLAZ_FLOWSOLVER_PARAM *_fp)
{
  _fp->nInnerElements.assign(M->z(),0);
  _fp->nBoundElements.assign(M->z(),0);
  for (size_t i=0; i<M->z(); ++i) {
    _fp->nInnerElements[i] = (M->d(i)==M->d()?   M->e(i) : 0);
    _fp->nBoundElements[i] = (M->d(i)==M->d()-1? M->e(i) : 0);
  }
  _fp->numDim = M->d();
  _fp->numUnk = M->d()+2;
  _fp->numNod = M->n();
  _fp->dtEul  = m_dt;

  _fp->time = 0.;
  _fp->iter = 0;
}
