
#include "mfactory.h"
#include "mfunction.h"
#include "t_manip.h"

using namespace m;


Register< mtransform,t_manip > mt_manip(12,"-tvsort","                variable sort by name (excluding coordinates)",
                                           "-tvkeep","[str:...]       ... to keep, removing the rest",
                                           "-tvrm",  "[str:...]       ... removal, by name",
                                           "-tvmv",  "[str=str:...]   ... moving, from name, to before name",
                                           "-tvren", "[str=str:...]   ... rename, from name, to name",
                                           "-tvadd", "[str[=fun]:...] ... addition, optionally set by function (uses only coordinate variables)",
                                           "-tvaxiz","                ... add new axisymmetric variables (around z axis, 3D only)",
                                           "-tzsort","                zone sort by name",
                                           "-tzkeep","[str:...]       ... to keep, removing the rest",
                                           "-tzrm",  "[str:...]       ... removal, by name",
                                           "-tzmv",  "[str=str:...]   ... moving, from name, to before name",
                                           "-tzren", "[str=str:...]   ... rename, from name, to name");


void t_manip::transform(GetPot& o, mmesh& m)
{
  using std::string;
  using std::vector;

  // get options key
  const string k = o[o.get_cursor()],
               v = (k=="-tvsort"||k=="-tvaxiz"||k=="-tzsort"? "" : o.get(o.inc_cursor(),""));

  // operations that apply in one shot
       if (k=="-tvsort") { vsort(m);   return; }
  else if (k=="-tzsort") { zsort(m);   return; }
  else if (k=="-tvkeep") { vkeep(m,v); return; }
  else if (k=="-tvaxiz") { vaxiz(m);   return; }
  else if (k=="-tzkeep") { zkeep(m,v); return; }

  // operations that apply multiple times
  const vector< std::pair< string,string > > vops = getoperands(v);
  for (vector< std::pair< string,string > >::const_iterator i=vops.begin(); i<vops.end(); ++i) {
         if (k=="-tvrm")  { vrm(m,getvindex(m,i->first)); }
    else if (k=="-tvmv")  { vmv(m,getvindex(m,i->first),getvindex(m,i->second)); }
    else if (k=="-tvren") { vren(m,getvindex(m,i->first),i->second); }
    else if (k=="-tvadd") { vadd(m,i->first,i->second,vector< string >(m.vn.begin(),m.vn.begin()+m.d())); }
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

  // create vector of (variable) names and sort it (excluding coordinates)
  vector< string > vn(m.vn);
  sort(vn.begin()+d,vn.end());

  // move (variable) name to new (variable) position
  for (unsigned i=d; i<v; ++i)
    vmv(m,getvindex(m,vn[i]),i);
}

void t_manip::vkeep(m::mmesh& m, const std::string& s)
{
  using std::string;
  using std::vector;

  // get list of (variable) names to keep
  vector< std::pair< string,string > > vops = getoperands(s);
  vector< string > keep;
  for (vector< std::pair< string,string > >::const_iterator i=vops.begin(); i!=vops.end(); ++i)
    keep.push_back(i->first);

  // remove variables not in that list (apply in reverse for performance)
  for (vector< string >::const_reverse_iterator i=m.vn.rbegin(); i!=m.vn.rend(); ++i)
    if (!std::count(keep.begin(),keep.end(),*i))
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

void t_manip::vadd(mmesh& m, const std::string& n, const std::string& f, const std::vector< std::string >& v)
{
  using std::string;
  using std::vector;

  // add a new variable (with zeros) in the end
  const unsigned Nnode = m.n();

  // find this variable index (if it doesn't exist set to the end)
  const unsigned v_idx = std::find(m.vn.begin(),m.vn.end(),n)!=m.vn.end()?
    std::distance(m.vn.begin(),std::find(m.vn.begin(),m.vn.end(),n)) :
    m.v();
  if (v_idx>=m.v()) {
    // variable name wasn't found, push it back
    m.vn.push_back(n);
    m.vv.push_back(vector< double >(Nnode,0.));
  }
  else {
    // variable name was found, it exists already
    m.vv[v_idx].assign(Nnode,0.);
  }
  if (!f.length())
    return;

  // set mfunction
  mfunction mf(f,v);

  // set variables for this function
  vector< unsigned > vindex(v.size(),0);
  for (unsigned i=0; i<v.size(); ++i)
    vindex[i] = getvindex(m,v[i]);

  // evaluate at each node
  vector< double >& veval = m.vv[v_idx];
  vector< double > c(v.size(),0.);
  for (unsigned i=0; i<Nnode; ++i) {
    for (unsigned j=0; j<v.size(); ++j)
      c[j] = m.vv[ vindex[j] ][i];
    veval[i] = mf.eval(&c[0]);
  }
}

void t_manip::vaxiz(mmesh& m)
{
  using namespace std;

  if (m.d()<3) {
    cerr << "error: axisymmetry variables only for 3D" << endl;
    throw 42;
  }

  // transform the coordinates
  const string c_rho  ("axi_"+m.vn[0]+"_rho"  ),  // name for coordinates rho
               c_theta("axi_"+m.vn[0]+"_theta");  // ... and theta
  vector< string > v(m.vn.begin(),m.vn.begin()+2);
  vadd(m,c_rho,   "sqrt("+v[0]+"*"+v[0]+"+"+v[1]+"*"+v[1]+")", v);
  vadd(m,c_theta, "atan2("+v[1]+","+v[0]+")",                  v);
  v.push_back(c_rho);
  v.push_back(c_theta);

  // transform the other vector fields (skipping coordinates)
  const vector< bool > visvector = m.vvectors();
  for (unsigned i=2; i<visvector.size(); ++i)
    if (visvector[i]) {
      vector< string > w(v);   // { Px,Py, Prho,Ptheta...
      w.push_back(m.vn[i+0]);  //   Vx,Vy }
      w.push_back(m.vn[i+1]);
      vadd(m,"axi_"+w[4]+"_rho",   "("+w[4]+"*"+w[0]+"+"+w[5]+"*"+w[1]+")/"+w[2], w);
      vadd(m,"axi_"+w[4]+"_theta", "("+w[5]+"*"+w[0]+"-"+w[4]+"*"+w[1]+")/"+w[2], w);
    }
}

void t_manip::zsort(mmesh& m)
{
  using namespace std;
  const unsigned z = m.z();
  if (!z)
    return;

  // create vector of (zone) names and sort it
  vector< string > zn(z);
  for (unsigned i=0; i<z; ++i)
    zn[i] = m.vz[i].n;
  sort(zn.begin(),zn.end());

  // move (zone) name to new (zone) position
  for (unsigned i=0; i<z; ++i)
    zmv(m,getzindex(m,zn[i]),i);
}

void t_manip::zkeep(m::mmesh& m, const std::string& s)
{
  using std::string;
  using std::vector;

  // get list of (zone) names to keep
  vector< std::pair< string,string > > vops = getoperands(s);
  vector< string > keep;
  for (vector< std::pair< string,string > >::const_iterator i=vops.begin(); i!=vops.end(); ++i)
    keep.push_back(i->first);

  // remove zones not in that list (apply in reverse for performance)
  for (vector< mzone >::const_reverse_iterator i=m.vz.rbegin(); i!=m.vz.rend(); ++i)
    if (!std::count(keep.begin(),keep.end(),i->n))
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


