
#include "mfactory.h"
#include "mfunction.h"
#include "t_manip.h"

using namespace m;


Register< mtransform,t_manip > mt_manip(11,"-tvsort","                variable sort by name (excluding coordinates)",
                                           "-tvkeep","[str:...]       ... to keep, removing the rest",
                                           "-tvrm",  "[str:...]       ... removal, by name",
                                           "-tvmv",  "[str=str:...]   ... moving, from name, to before name",
                                           "-tvren", "[str=str:...]   ... rename, from name, to name",
                                           "-tvadd", "[str[=fun]:...] ... addition, optionally set by function",
                                           "-tzsort","                zone sort by name",
                                           "-tzkeep","[str:...]       ... to keep, removing the rest",
                                           "-tzrm",  "[str:...]       ... removal, by name",
                                           "-tzmv",  "[str=str:...]   ... moving, from name, to before name",
                                           "-tzren", "[str=str:...]   ... rename, from name, to name");


void t_manip::transform(GetPot& o, mmesh& m)
{
  using std::string;

  // get options key
  const string k = o[o.get_cursor()],
               v = (k=="-tvsort"||k=="-tzsort"? "" : o.get(o.inc_cursor(),""));

  // operations that apply in one shot
       if (k=="-tvsort") { vsort(m);   return; }
  else if (k=="-tzsort") { zsort(m);   return; }
  else if (k=="-tvkeep") { vkeep(m,v); return; }
  else if (k=="-tzkeep") { zkeep(m,v); return; }

  // operations that apply multiple times
  const std::vector< std::pair< string,string > > vops = getoperands(v);
  for (std::vector< std::pair< string,string > >::const_iterator i=vops.begin(); i<vops.end(); ++i) {
         if (k=="-tvrm")  { vrm(m,getvindex(m,i->first)); }
    else if (k=="-tvmv")  { vmv(m,getvindex(m,i->first),getvindex(m,i->second)); }
    else if (k=="-tvren") { vren(m,getvindex(m,i->first),i->second); }
    else if (k=="-tvadd") { vadd(m,i->first,i->second); }
    else if (k=="-tzrm")  { zrm(m,getzindex(m,i->first)); }
    else if (k=="-tzmv")  { zmv(m,getzindex(m,i->first),getzindex(m,i->second)); }
    else if (k=="-tzren") { zren(m,getzindex(m,i->first),i->second); }
  }
}

void t_manip::vsort(mmesh& m)
{
  using namespace std;
  const unsigned d = m.d(),
                 v = m.v();
  if (v<=d)
    return;

  // create vector of names / (current) indices
  vector< pair< string,unsigned > > vn(v);
  for (unsigned i=d; i<v; ++i)
    vn[i] = make_pair(m.vn[i],i);

  // sort by name and move (variable) name to new (variable) position
  sort(vn.begin()+d,vn.end());
  for (unsigned i=d; i<v; ++i)
    vmv(m,getvindex(m,vn[i].first),i);
}

void t_manip::vkeep(m::mmesh& m, const std::string& s)
{
  using std::string;
  using std::vector;

  // get list of (variable) names to keep
  vector< std::pair< string,string > > vops = getoperands(s);
  vector< string > vkeep;
  for (vector< std::pair< string,string > >::const_iterator i=vops.begin(); i!=vops.end(); ++i)
    vkeep.push_back(i->first);

  // remove variables not in that list (apply in reverse for performance)
  for (vector< string >::const_reverse_iterator i=m.vn.rbegin(); i!=m.vn.rend(); ++i)
    if (!std::count(vkeep.begin(),vkeep.end(),*i))
      vrm(m,getvindex(m,*i));
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

void t_manip::vren(mmesh& m, const unsigned i, const std::string& n)
{
  m.vn[i] = n;
}

void t_manip::vadd(mmesh& m, const std::string& n, const std::string& f)
{
  // add a new variable (with zeros) in the end
  const unsigned Ndim  = m.d();
  const unsigned Nnode = m.n();
  m.vn.push_back(n);
  m.vv.push_back(std::vector< double >(Nnode,0.));
  if (!f.length())
    return;

  // set mfunction (uses only the coordinates variables, for your own safety)
  const std::vector< std::string > vnames(m.vn.begin(),m.vn.begin()+Ndim);
  mfunction mf(f,vnames);

  // evaluate at each node
  std::vector< double >& v = m.vv.back();
  std::vector< double > c(Ndim,0.);
  for (unsigned n=0; n<Nnode; ++n) {
    for (unsigned d=0; d<Ndim; ++d)
      c[d] = m.vv[d][n];
    v[n] = mf.eval(&c[0]);
  }
}

void t_manip::zsort(mmesh& m)
{
  using namespace std;
  const unsigned z = m.z();
  if (!z)
    return;

  // create vector of names / (current) indices
  vector< pair< string,unsigned > > zn(z);
  for (unsigned i=0; i<z; ++i)
    zn[i] = make_pair(m.vz[i].n,i);

  // sort by name and move (zone) name to new (zone) position
  sort(zn.begin(),zn.end());
  for (unsigned i=0; i<z; ++i)
    zmv(m,getzindex(m,zn[i].first),i);
}

void t_manip::zkeep(m::mmesh& m, const std::string& s)
{
  using std::string;
  using std::vector;

  // get list of (zone) names to keep
  vector< std::pair< string,string > > vops = getoperands(s);
  vector< string > zkeep;
  for (vector< std::pair< string,string > >::const_iterator i=vops.begin(); i!=vops.end(); ++i)
    zkeep.push_back(i->first);

  // remove zones not in that list (apply in reverse for performance)
  for (vector< mzone >::const_reverse_iterator i=m.vz.rbegin(); i!=m.vz.rend(); ++i)
    if (!std::count(zkeep.begin(),zkeep.end(),i->n))
      zrm(m,getzindex(m,i->n));
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

void t_manip::zren(mmesh& m, const unsigned i, const std::string& n)
{
  m.vz[i].n = n;
}

unsigned t_manip::getvindex(const mmesh& m, const std::string& n)
{
  for (unsigned i=0; i<m.vn.size(); ++i)
    if (n==m.vn[i])
      return i;
  std::cerr << "error: variable name not found: \"" << n << "\"" << std::endl;
  throw 42;
  return 0;
}

unsigned t_manip::getzindex(const mmesh& m, const std::string& n)
{
  for (unsigned i=0; i<m.vz.size(); ++i)
    if (n==m.vz[i].n)
      return i;
  std::cerr << "error: zone name not found: \"" << n << "\"" << std::endl;
  throw 42;
  return 0;
}

std::vector< std::pair< std::string,std::string > >
  t_manip::getoperands(const std::string& s)
{
  using std::string;

  // split string find a '=' then a ':'
  std::vector< std::pair< string,string > > r;
  string::size_type p1 = 0,
                         p2 = 0;
  while (p2!=string::npos) {
    std::pair< string,string > p("","");

    p2 = std::min(s.find(":",p1),s.find("=",p1));
    p.first = s.substr(p1,(p2==string::npos? p2:p2-p1));
    if (s.find("=",p1)<s.find(":",p1)) {
      p1 = p2+1;
      p2 = s.find(":",p1);

      p.second = s.substr(p1,(p2==string::npos? p2:p2-p1));
      /* older version, maybe works better
        string s2 = s.substr(p1,(p2==string::npos? p2:p2-p1));
        istringstream ss(s2);
        ss >> p.second;
      */
    }

    p1 = p2+1;
    r.push_back(p);
  }
  return r;
}


