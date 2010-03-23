
#include <sstream>
#include "mfactory.h"
#include "mfunction.h"
#include "t_manip.h"

using namespace std;
using namespace m;


Register< mtransform,t_manip > mt_manip(7,"-tvrm", "[str:...] variable removal, by name",
                                          "-tvmv", "[str=str:...]   variable moving, from name, to before name",
                                          "-tvren","[str=str:...]   variable rename, from name, to name",
                                          "-tvadd","[str[=fun]:...] variable addition, optionally set by function",
                                          "-tzrm", "[str:...] zone removal, by name",
                                          "-tzmv", "[str=str:...]   zone moving, from name, to before name",
                                          "-tzren","[str=str:...]   zone rename, from name, to name");


void t_manip::transform(GetPot& o, mmesh& m)
{
  // get key and value
  const string k = o[o.get_cursor()];
  const vector< pair< string,string > > vops = getoperations(o.get(o.inc_cursor(),""));
  for (vector< pair< string,string > >::const_iterator i=vops.begin(); i<vops.end(); ++i) {
         if (k=="-tvrm")  { vrm(m,getvindex(m,i->first)); }
    else if (k=="-tvmv")  { vmv(m,getvindex(m,i->first),getvindex(m,i->second)); }
    else if (k=="-tvren") { vren(m,getvindex(m,i->first),i->second); }
    else if (k=="-tvadd") { vadd(m,i->first,i->second); }
    else if (k=="-tzrm")  { zrm(m,getzindex(m,i->first)); }
    else if (k=="-tzmv")  { zmv(m,getzindex(m,i->first),getzindex(m,i->second)); }
    else if (k=="-tzren") { zren(m,getzindex(m,i->first),i->second); }
  }
}

void t_manip::vrm(mmesh& m, const unsigned i)
{
  m.vn.erase(m.vn.begin()+i);
  m.vv.erase(m.vv.begin()+i);
}

void t_manip::vmv(mmesh& m, const unsigned i, const unsigned j)
{
  m.vn.insert(m.vn.begin()+j,*(m.vn.begin()+i));
  m.vv.insert(m.vv.begin()+j,*(m.vv.begin()+i));
  vrm(m,i<j? i:i+1);
}

void t_manip::vren(mmesh& m, const unsigned i, const string& n)
{
  m.vn[i] = n;
}

void t_manip::vadd(mmesh& m, const string& n, const string& f)
{
  // add a new variable (with zeros) in the end
  const unsigned Ndim  = m.d();
  const unsigned Nnode = m.n();
  m.vn.push_back(n);
  m.vv.push_back(vector< double >(Nnode,0.));
  if (!f.length())
    return;

  // set mfunction (uses only the coordinates variables, for your own safety)
  const vector< string > vnames(m.vn.begin(),m.vn.begin()+Ndim);
  mfunction mf(f,vnames);

  // evaluate at each node
  vector< double >& v = m.vv.back();
  vector< double > c(Ndim,0.);
  for (unsigned n=0; n<Nnode; ++n) {
    for (unsigned d=0; d<Ndim; ++d)
      c[d] = m.vv[d][n];
    v[n] = mf.eval(&c[0]);
  }
}

void t_manip::zrm(mmesh& m, const unsigned i)
{
  m.vz.erase(m.vz.begin()+i);
}

void t_manip::zmv(mmesh& m, const unsigned i, const unsigned j)
{
  m.vz.insert(m.vz.begin()+j,*(m.vz.begin()+i));
  zrm(m,i<j? i:i+1);
}

void t_manip::zren(mmesh& m, const unsigned i, const string& n)
{
  m.vz[i].n = n;
}

unsigned t_manip::getvindex(const m::mmesh& m, const string& n)
{
  for (unsigned i=0; i<m.vn.size(); ++i)
    if (n==m.vn[i])
      return i;
  cerr << "error: variable name not found: \"" << n << "\"" << endl;
  throw 42;
  return 0;
}

unsigned t_manip::getzindex(const m::mmesh& m, const string& n)
{
  for (unsigned i=0; i<m.vz.size(); ++i)
    if (n==m.vz[i].n)
      return i;
  cerr << "error: zone name not found: \"" << n << "\"" << endl;
  throw 42;
  return 0;
}

vector< pair< string,string > > t_manip::getoperations(const string& s)
{
  // split string find a '=' then a ':'
  vector< pair< string,string > > r;
  string::size_type p1 = 0,
                    p2 = 0;
  while (p2!=string::npos) {
    pair< string,string > p("","");

    p2 = min(s.find(":",p1),s.find("=",p1));
    p.first = s.substr(p1,(p2==string::npos? p2:p2-p1));
    if (s.find("=",p1)<s.find(":",p1)) {
      p1 = p2+1;
      p2 = s.find(":",p1);
      string s2 = s.substr(p1,(p2==string::npos? p2:p2-p1));
      istringstream ss(s2);
      ss >> p.second;
    }

    p1 = p2+1;
    r.push_back(p);
  }
  return r;
}


