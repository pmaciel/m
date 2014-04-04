
#include <sstream>
#include "mfactory.h"
#include "f_xpl.h"

using namespace std;
using namespace m;


Register< mfinput,f_xpl > mf_xpl(".xpl","xplot input format");


void f_xpl::read(GetPot& o, mmesh& m)
{
  const string fn(o.get(o.inc_cursor(),""));
  ifstream f(fn.c_str());
  if (!f) {
    cerr << "error accessing file: \"" << fn << "\"" << endl;
    throw 42;
  }
  string s;      // auxiliary string
  unsigned dum;  // ... unsigned integer

  // skip first line
  getline(f,s);

  // allocate memory for node coordinates
  unsigned Nnode, Nelem, Nbelem, Nvar;
  {
    getline(f,s);
    istringstream ss(s);
    ss >> Nelem >> Nnode >> Nbelem >> Nvar;
    m.vv.assign(2+Nvar,vector< double >(Nnode,0.));
    m.vn.resize(2+Nvar);
    m.vn[0] = "x";
    m.vn[1] = "y";
    for (unsigned i=0; i<Nvar; ++i)
      m.vn[i+2] = "Var" + string(1,char('A'+i)) + '_';
  }

  // set first (inner) zone
  melem default_e;
  mzone default_z;
  default_e.n.assign(3,0);
  default_z.n = "InnerCells";
  default_z.e2n.assign(Nelem,default_e);
  default_z.t = FETRIANGLE;
  m.vz.assign(1,default_z);
  m.vz[0].n = "InnerCells";
  m.vz[0].e2n.assign(Nelem,default_e);
  m.vz[0].t = FETRIANGLE;

  // read elements connectivity
  for (unsigned i=0; i<Nelem; ++i) {
    getline(f,s);
    istringstream ss(s);
    ss >> dum;
    vector< unsigned >& en = m.vz[0].e2n[i].n;
    for (unsigned j=0; j<en.size(); ++j) {
      ss >> en[j];
      --en[j];
    }
  }

  // skip number of nodes (again?) and "dummy free-stream values"
  getline(f,s);
  getline(f,s);
  for (unsigned i=0; i<Nnode; ++i) {
    getline(f,s);
    istringstream ss(s);
    for (unsigned j=0; j<2+Nvar; ++j)
      ss >> m.vv[j][i];
  }

  // skip number of boundaries and boundary elements (again?)
  getline(f,s);
  getline(f,s);

  // read boundaries into a matrix [col][line] to convert later
  vector< vector< unsigned > > b(3,vector< unsigned >(Nbelem,0));
  for (unsigned i=0; i<Nbelem; ++i) {
    getline(f,s);
    istringstream ss(s);
    ss >> b[0][i] >> b[1][i] >> b[2][i];
  }

  // resize boundary zones
  const unsigned Nbzones = *max_element(b[2].begin(),b[2].end());
  default_e.n.assign(2,0);
  default_z.n = "Boundary";  // later corrected with a number
  default_z.e2n.resize(0);   // later resized
  default_z.t = FELINESEG;
  m.vz.resize(1+Nbzones,default_z);
  for (unsigned i=0; i<Nbzones; ++i) {
    const unsigned Nfaces = count_if(b[2].begin(),b[2].end(),
      bind2nd(equal_to< unsigned >(),i+1) );
    m.vz[i+1].e2n.reserve(Nfaces);
    m.vz[i+1].n += string(1,char('1'+i));
  }

  // convert boundary zones matrix
  for (unsigned i=0; i<Nbelem; ++i) {
    default_e.n[0] = b[0][i]-1;
    default_e.n[1] = b[1][i]-1;
    m.vz[ b[2][i] ].e2n.push_back(default_e);
  }

  // skip next line (what is it, '"'?)
  getline(f,s);
  cout << '"' << s << '"' << endl;

  // read a strange field, i think it should be stored somewhere else
  // (in case of reading it successfully, store it to the last variable)
  while (true) {
    vector< double > v;
    v.reserve(Nnode);
    double value;  // or your money back
    while (f >> value)
      v.push_back(value);
    if (v.size()==Nnode) {
      m.vv.push_back(v);
      m.vn.push_back("Var" + string(1,char('A'+(Nvar++))) + '_');
    }
    else
      break;
  }

  // close file
  f.close();
}

