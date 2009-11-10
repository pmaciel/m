
#include <sstream>
#include "mfactory.h"
#include "f_cfmesh.h"

using namespace std;
using namespace m;


Register< mfinput,f_cfmesh > mf_cfmesh1(2,".CFmesh","CFmesh input format",
                                          ".CFmeshFE","(PlatingMaster extension variant)");
Register< mfoutput,f_cfmesh > mf_cfmesh2(".CFmesh","CFmesh output format");


void f_cfmesh::read(GetPot& o, mmesh& m)
{
  const string fn(o.get(o.inc_cursor(),""));
  ifstream f(fn.c_str());
  if (!f) {
    cerr << "error accessing file: \"" << fn << "\"" << endl;
    throw 42;
  }

  // auxiliary variables
  istringstream ss;
  string s, TRS_NAME;
  unsigned NB_DIM=0, NB_EQ=0,                      // general setup
           NB_NODES=0, NB_STATES=0, LIST_STATE=0,  // nodal setup
           NB_ELEM=0, NB_NODES_PER_TYPE=0,         // volume connectivity setup
           NB_TRSs=0, NB_GEOM_ENTS=0;              // boundary connectivities setup

  while (!f.eof()) {
    getline(f,s);  ss.clear();  ss.str(s);  ss >> s;

    // general setup: version, exit dimensions and variables number
    if      (s=="!COOLFLUID_VERSION")     {  ss >> s;  cout << "info: coolfluid version: "    << s << endl;  }
    else if (s=="!COOLFLUID_SVNVERSION")  {  ss >> s;  cout << "info: coolfluid svnversion: " << s << endl;  }
    else if (s=="!CFMESH_FORMAT_VERSION") {  ss >> s;  cout << "info: format version: " << s << endl;  }
    else if (s=="!END")     {  break;  }
    else if (s=="!NB_DIM")  {  ss >> NB_DIM;  if (NB_DIM && NB_EQ)  setvarnames(m.vn,NB_DIM,NB_EQ);  }
    else if (s=="!NB_EQ")   {  ss >> NB_EQ;   if (NB_DIM && NB_EQ)  setvarnames(m.vn,NB_DIM,NB_EQ);  }

    // nodal setup: nodes (coordinates) and states (variable values)
    else if (s=="!NB_NODES")    {  ss >> NB_NODES;   }
    else if (s=="!NB_STATES")   {  ss >> NB_STATES;  }
    else if (s=="!LIST_NODE")   {
      if (NB_NODES)  readlnodes(f,m,NB_DIM,NB_EQ,NB_NODES);
    }
    else if (s=="!LIST_STATE")  {
      ss >> LIST_STATE;
      if      (NB_STATES && LIST_STATE)   readlstates1(f,m,NB_DIM,NB_EQ,NB_STATES);
      else if (NB_STATES && !LIST_STATE)  readlstates0(f,m,NB_DIM,NB_EQ,NB_STATES);
    }

    // volume connectivity setup
    else if (s=="!NB_ELEM")             {}
    else if (s=="!GEOM_POLYORDER")      {}
    else if (s=="!SOL_POLYORDER")       {}
    else if (s=="!NB_NODES_PER_TYPE")   {  ss >> NB_NODES_PER_TYPE;   }
    else if (s=="!ELEM_TYPES")          {}
    else if (s=="!NB_STATES_PER_TYPE")  {}
    else if (s=="!NB_ELEM_TYPES")       {}
    else if (s=="!NB_ELEM_PER_TYPE")    {
      unsigned Nelem = 0;
      while (ss >> NB_ELEM)
        Nelem += NB_ELEM;
      NB_ELEM = Nelem;
    }
    else if (s=="!LIST_ELEM") {

      // set zone header
      m.vz.resize(1);
      mzone& z = m.vz.back();
      z.n = "InnerCells";
      z.t = (NB_DIM==2? FETRIANGLE : (NB_DIM==3? FETETRAHEDRON : ORDERED ));

      // set zone connectivity list
      if (NB_ELEM && NB_NODES_PER_TYPE) {
        melem e;
        e.n.assign(NB_NODES_PER_TYPE,0);
        z.e2n.assign(NB_ELEM,e);
        readlelems(f,z.e2n,NB_NODES_PER_TYPE,NB_ELEM);
      }

    }

    // boundary connectivities setup
    else if (s=="!NB_TRSs")        {  ss >> NB_TRSs;   }
    else if (s=="!TRS_NAME")       {  ss >> TRS_NAME;  }
    else if (s=="!NB_TRs")         {}
    else if (s=="!GEOM_TYPE")      {}
    else if (s=="!NB_GEOM_ENTS")   {
      unsigned Nelem = 0;
      while (ss >> NB_GEOM_ENTS)
        Nelem += NB_GEOM_ENTS;
      NB_GEOM_ENTS = Nelem;
    }
    else if (s=="!LIST_GEOM_ENT") {

      // set zone header
      m.vz.reserve(1+NB_TRSs);
      m.vz.push_back(mzone());
      mzone& z = m.vz.back();
      z.n = TRS_NAME;
      z.t = (NB_DIM==2? FELINESEG : (NB_DIM==3? FETRIANGLE : ORDERED ));

      // set zone connectivity list
      if (NB_GEOM_ENTS) {
        melem e;
        e.n.assign(NB_DIM,0);
        z.e2n.assign(NB_GEOM_ENTS,e);
        readlelems(f,z.e2n,0,NB_GEOM_ENTS);
      }

    }
    else {
      cout << "expecting key, found \"" << ss.str() << "\"" << endl;
    }
  }

#if 0
  // flip 3d boundary elements
  if (NB_DIM==3) {
    for (unsigned t=1; t<m.vz.size(); ++t)
      for (unsigned c=0; c<m.e(t); ++c)
        swap(m.vz[t].e2n[c].n[0],m.vz[t].e2n[c].n[1]);
  }
#endif
}


void f_cfmesh::write(GetPot& o, const mmesh& m)
{
  const string fn(o.get(o.inc_cursor(),""));
  ofstream f(fn.c_str());
  if (!f) {
    cerr << "error accessing file: \"" << fn << "\"" << endl;
    throw 42;
  }


  // reorder zones, setting the first one as the volume

  const unsigned dim = m.d();
  vector< const mzone* > z;
  for (unsigned i=0; i<m.vz.size(); ++i) {
    if (m.d(i)==dim && (m.vz[i].t==FETRIANGLE || m.vz[i].t==FETETRAHEDRON)) {
      z.push_back(&m.vz[i]);
      break;  // only one
    }
  }
  if (z.size()!=1) {
    cerr << "couldn't find a volume zone!" << endl;
    throw 42;
  }
  for (unsigned i=0; i<m.vz.size(); ++i) {
    if (m.d(i)==dim-1)
      z.push_back(&m.vz[i]);
  }


  // header
  f << "!COOLFLUID_VERSION 2.0.0" << endl
    << "!COOLFLUID_SVNVERSION 8703M" << endl
    << "!CFMESH_FORMAT_VERSION 1.3" << endl;
  f << "!NB_DIM " << dim << endl
    << "!NB_EQ " << max< unsigned >(m.v()-dim,1) << endl
    << "!NB_NODES " << m.n() << " 0" << endl
    << "!NB_STATES " << m.n() << " 0" << endl;


  // elements connectivity, volume (0-based)
  f << "!NB_ELEM " << m.e() << endl
    << "!NB_ELEM_TYPES 1" << endl
    << "!GEOM_POLYORDER 1" << endl
    << "!SOL_POLYORDER 1" << endl
    << "!ELEM_TYPES " << (dim>2? "Tetra":"Triag") << endl
    << "!NB_ELEM_PER_TYPE " << z[0]->e2n.size() << endl
    << "!NB_NODES_PER_TYPE " << dim+1 << endl
    << "!NB_STATES_PER_TYPE " << dim+1 << endl;
  f << "!LIST_ELEM" << endl;
  for (unsigned c=0; c<z[0]->e2n.size(); ++c) {
    const vector< unsigned >& en = z[0]->e2n[c].n;
    for (unsigned i=0; i<en.size(); ++i)
      f << " " << en[i];
    for (unsigned i=0; i<en.size(); ++i)
      f << " " << en[i];
    f << endl;
  }


  // elements connectivity, boundaries (0-based)
  f << "!NB_TRSs " << z.size()-1 << endl;
  for (unsigned t=1; t<z.size(); ++t) {
    f << "!TRS_NAME " << z[t]->n << endl
      << "!NB_TRs 1" << endl
      << "!NB_GEOM_ENTS " << z[t]->e2n.size() << endl
      << "!GEOM_TYPE Face" << endl
      << "!LIST_GEOM_ENT" << endl;
    for (unsigned c=0; c<z[t]->e2n.size(); ++c) {
      const vector< unsigned >& en = z[t]->e2n[c].n;
      f << en.size() << ' ' << en.size();
      for (unsigned i=0; i<en.size(); ++i)
        f << ' ' << en[i];
      for (unsigned i=0; i<en.size(); ++i)
        f << ' ' << en[i];
      f << endl;
    }
  }


  // coordinates
  f << "!LIST_NODE" << endl;
  for (unsigned n=0; n<m.n(); ++n) {
    for (unsigned d=0; d<dim; ++d)
      f << ' ' << m.vv[d][n];
    f << endl;
  }


  // states
  f << "!LIST_STATE " << (m.v()>dim? "1":"0") << endl;
  if (m.v()>dim)
    for (unsigned n=0; n<m.n(); ++n) {
      for (unsigned d=dim; d<m.v(); ++d)
        f << " " << m.vv[d][n];
      f << endl;
    }


  f << "!END" << endl;
  f.close();
}


void f_cfmesh::setvarnames(vector< string >& vn, unsigned Ndim, unsigned Neqs)
{
  vn.resize(Ndim+Neqs);
  for (unsigned d=0; d<Ndim; ++d) {
    vn[d] = string(1,char('x'+d));
  }
  for (unsigned d=0; d<Neqs; ++d) {
    ostringstream ss;  ss << "Var" << char('A'+d) << '_';  vn[Ndim+d] = ss.str();
  }
}


void f_cfmesh::readlnodes(ifstream& f, mmesh& m, unsigned Ndim, unsigned Neqs, unsigned N)
{
  if (!m.vv.size())
    m.vv.assign(Ndim+Neqs,vector< double >(N,0.));

  string s;
  istringstream ss;

  for (unsigned n=0; n<N; ++n) {
    getline(f,s);  ss.clear();  ss.str(s);
    for (unsigned d=0; d<Ndim; ++d)
      ss >> m.vv[d][n];
  }
}


void f_cfmesh::readlstates1(ifstream& f, mmesh& m, unsigned Ndim, unsigned Neqs, unsigned N)
{
  if (!m.vv.size())
    m.vv.assign(Ndim+Neqs,vector< double >(N,0.));

  istringstream ss;
  string s;
  for (unsigned n=0; n<N; ++n) {
    getline(f,s);  ss.clear();  ss.str(s);
    for (unsigned d=Ndim; d<Ndim+Neqs; ++d)
      ss >> m.vv[d][n];
  }
}


void f_cfmesh::readlstates0(ifstream& f, mmesh& m, unsigned Ndim, unsigned Neqs, unsigned N)
{
  if (!m.vv.size())
    m.vv.assign(Ndim,vector< double >(N,0.));

  vector< vector< double > > vvkeep(Ndim);
  vector< string >           vnkeep(Ndim);
  for (unsigned d=0; d<Ndim; ++d) {
    m.vv[d].swap(vvkeep[d]);
    m.vn[d].swap(vnkeep[d]);
  }
  m.vv.swap(vvkeep);
  m.vn.swap(vnkeep);
}


void f_cfmesh::readlelems(ifstream& f, vector< melem >& e2n, unsigned Nnodes, unsigned N)
{
  istringstream ss;
  string s;
  unsigned nnodes  = Nnodes,
           nstates = Nnodes;
  for (vector< melem >::iterator c=e2n.begin(); c!=e2n.end(); ++c) {
    getline(f,s);  ss.clear();  ss.str(s);
    if (!Nnodes) {
      ss >> nnodes >> nstates;
      c->n.resize(nnodes);
    }
    for (unsigned i=0; i<nnodes; ++i)
      ss >> c->n[i];
  }
}

