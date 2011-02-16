
//include <algorithm>  // for transform
//include <cctype>     // for toupper
#include <iostream>
#include <sstream>
#include "mfactory.h"
#include "f_plt.h"

using namespace std;
using namespace m;


Register< mfinput,f_plt > mf_plt1(3,".plt","Tecplot ASCII input format",
                                    "","--plt-pm2d: input is PlatingMaster 2D format (default: no)",
                                    "","--plt-pm3d: input is PlatingMaster 3D format (default: no)");
Register< mfoutput,f_plt > mf_plt2(3,".plt","Tecplot ASCII output format",
                                     "","--plt-block: output in point format (default: no)",
                                     "","--plt-solutiontime [real]: solution time (default: 0.)");


void f_plt::read(GetPot& o, mmesh& m)
{
  const string fn(o.get(o.inc_cursor(),""));
  ifstream f(fn.c_str());
  if (!f) {
    cerr << "error accessing file: \"" << fn << "\"" << endl;
    throw 42;
  }

  if      (o.search("--plt-pm2d")) { in_platingmaster_2d(f,m); }
  else if (o.search("--plt-pm3d")) { in_platingmaster_3d(f,m); }
  else {

    string s;
    while (getline(f,s)) {
      const string key(upper(splitstring(s)[0]));
      if (key.find("VARIABLES")==0) {

        m.vn = getVariables(s);

      }
      else if (key=="ZONE") {

        m.vz.push_back(mzone());
        mzone& z = m.vz.back();
        TecZone tz = getZoneHeader(s,z.n,z.t);
        if (!tz.isshared)
          readZoneNodeValues(f,m.vv,tz.i,m.v(),tz.isblock);
        if (tz.nenodes)
          readZoneConnectivity(f,z.e2n,tz.j,tz.nenodes,z.t);

      }
    }

  }

  f.close();
}


void f_plt::write(GetPot& o, const mmesh& m)
{
  const string fn(o.get(o.inc_cursor(),""));
  ofstream f(fn.c_str());
  if (!f) {
    cerr << "error accessing file: \"" << fn << "\"" << endl;
    throw 42;
  }
  f.precision(15);
  const bool   isblock      = o.search("--plt-block");
  const double solutiontime = o.follow(0.,"--plt-solutiontime");


  // header, then zone header, nodal and elements information
  f << setVariables(m.vn) << endl;
  for (unsigned t=0; t<m.vz.size(); ++t) {
    const mzone& z = m.vz[t];
    switch (z.t) {

      // these are the naturally supported element types
      case ORDERED:
      case FELINESEG:
      case FETRIANGLE:
      case FEQUADRILATERAL:
      case FETETRAHEDRON:
      case FEBRICK:
      case FEPOLYGON:
      case FEPOLYHEDRON: {

        TecZone tz(z.n,                                // t
                   z.e2n.size()? z.e2n[0].n.size():0,  // nenodes
                   m.n(),                              // i
                   z.e2n.size(),                       // j
                   z.t,                                // et
                   isblock,                            // isblock
                   t>0);                               // isshared
        f << setZoneHeader(tz,m.v(),solutiontime) << endl;
        writeZoneNodeValues(f,m.vv,isblock,t>0? 0:-1);  // (isblock is an option)
        writeZoneConnectivity(f,z.e2n,z.t);

      } break;

      // these are represented by a FEBRICK, with coalesced nodes
      case PRISM3:
      case PYRAMID4: {

        melem e;
        e.n.assign(8,0);
        vector< melem > e2n(z.e2n.size(),e);
        for (unsigned c=0; c<z.e2n.size(); ++c) {
          const vector< unsigned >& eng = z.e2n[c].n;  // element nodes, original
                vector< unsigned >& en  =   e2n[c].n;  // ...,           modified
          if (z.t==PRISM3) {
            en[0] = eng[0];
            en[1] = eng[1];
            en[2] = eng[2];
            en[3] = eng[2];
            en[4] = eng[3];
            en[5] = eng[4];
            en[6] = eng[5];
            en[7] = eng[5];
          }
          else if (z.t==PYRAMID4) {
            en[0] = eng[0];
            en[1] = eng[1];
            en[2] = eng[2];
            en[3] = eng[3];
            en[4] = eng[4];
            en[5] = eng[4];
            en[6] = eng[4];
            en[7] = eng[4];
          }
        }
        mtype et = FEBRICK;
        TecZone tz(z.n,         // t
                   z.d()+1,     // nenodes
                   m.n(),       // i
                   e2n.size(),  // j
                   et,          // et
                   false,       // isblock
                   t>0);        // isshared
        f << setZoneHeader(tz,m.v(),solutiontime) << endl;
        writeZoneNodeValues(f,m.vv,isblock,t>0? 0:-1);
        writeZoneConnectivity(f,e2n,et);

      } break;
      default: break;

    }
  }
  f.close();
}


// PlatingMaster specific functions ------------------------------------------


void f_plt::in_platingmaster_2d(ifstream& f, mmesh& m)
{
  string s;
  bool discardx = false,
       discardy = false,
       discardz = false;

  while (getline(f,s)) {
    const string key(upper(splitstring(s)[0]));
    if (key.find("VARIABLES")==0) {

      m.vn = getVariables(s);
      if (m.d()!=3) {
        cerr << "error: PlatingMaster unexpected number of dimensions!" << endl;
        throw 42;
      }

    }
    else if (key=="ZONE") {

      m.vz.push_back(mzone());
      mzone& z = m.vz.back();
      TecZone tz = getZoneHeader(s,z.n,z.t);
      while (z.n.find(" ")!=string::npos)
        z.n.replace(z.n.find(" "),1,"_");

      if (m.vz.size()==1) {

        // read a FACE grid
        z.t = FETRIANGLE;
        readZoneNodeValues(f, m.vv, tz.i,3,false); // there is always a 3rd coordinate, to be discarded
        readZoneConnectivity(f,z.e2n,tz.j,tz.nenodes,z.t);

        // guess coordinate to discard and make sure z is empty
        vector< double >& vx = m.vv[0];
        vector< double >& vy = m.vv[1];
        vector< double >& vz = m.vv[2];
        const double dx = *max_element(vx.begin(),vx.end()) - *min_element(vx.begin(),vx.end());
        const double dy = *max_element(vy.begin(),vy.end()) - *min_element(vy.begin(),vy.end());
        const double dz = *max_element(vz.begin(),vz.end()) - *min_element(vz.begin(),vz.end());
        discardx = (dx<dy && dx<dz);
        discardy = (dy<dx && dy<dz);
        discardz = (dz<dx && dz<dy);
        if      (discardx)  { m.vv[0].swap(m.vv[2]); }
        else if (discardy)  { m.vv[1].swap(m.vv[2]); }
        else if (discardz)  {}
        else {
          cerr << "error: couldn't distinguish a variable to discard" << endl;
          throw 42;
        }
        m.vv.pop_back();
        m.vn.pop_back();

        // swap triangle nodes because PM convention have a "positive" z
        for (unsigned c=0; c<m.e() && z.e2n[c].n.size()==3; ++c) {
          vector< unsigned >& en = z.e2n[c].n;
          const double v1x = m.vv[0][en[1]]-m.vv[0][en[0]],
                       v1y = m.vv[1][en[1]]-m.vv[1][en[0]],
                       v2x = m.vv[0][en[2]]-m.vv[0][en[0]],
                       v2y = m.vv[1][en[2]]-m.vv[1][en[0]];
          if (v1x*v2y-v1y*v2x<0.)
            swap(en[0],en[1]);
        }

      }
      else {

        // read an EDGE grid
        z.t = FELINESEG;
        vector< vector< double > > vv;
        readZoneNodeValues(f, vv, tz.i,3,false); // there is always a 3rd coordinate, to be discarded
        if      (discardx)  { vv[0].swap(vv[2]); }
        else if (discardy)  { vv[1].swap(vv[2]); }
        else if (discardz)  {}
        else {
          cerr << "error: couldn't distinguish a variable to discard" << endl;
          throw 42;
        }
        vv.pop_back();

        cout << "info: reconstruct boundary zone " << m.vz.size()-1 << "..." << endl;
        vector< unsigned > vi;
        getPMBoundaryZonePointNodeIndices(m,vi,vv[0],vv[1]);
        getPMBoundaryZoneElementsFromNodeIndices(m,vi,vv[0],vv[1], z.e2n);

      }

    }
  }
}


void f_plt::in_platingmaster_3d(ifstream& f, mmesh& m)
{
  cerr << "error: PlatingMaster 3 dimensions not implemented!" << endl;
  throw 42;
}


void f_plt::getPMBoundaryZonePointNodeIndices(const mmesh& m, vector< unsigned >& vi, const vector< double >& vx, const vector< double >& vy)
{
  // compare (square of) distances
  vi.resize(vx.size());
  for (unsigned j=0; j<vx.size(); ++j) {
    unsigned ni = 42;   // closest node index
    double nd = 1.e99;  // ... distance
    for (unsigned n=0; n<m.n(); ++n) {
      const double d = (vx[j]-m.vv[0][n])*(vx[j]-m.vv[0][n]) +
                       (vy[j]-m.vv[1][n])*(vy[j]-m.vv[1][n]);
      if (d<nd) {
        nd=d;
        ni=n;
      }
    }
    vi[j] = ni;
  }
}


void f_plt::getPMBoundaryZoneElementsFromNodeIndices(const mmesh& m, const vector< unsigned >& vi, const vector< double >& vx, const vector< double >& vy, vector< melem >& ve)
{
  const vector< melem >& alle = m.vz[0].e2n;  // all elements

  // create list of elements with any of these nodes
  vector< unsigned > vve;  // vector of volume elements
  for (unsigned i=0; i<vi.size(); ++i)
    for (unsigned c=0; c<alle.size(); ++c)
      if (count(alle[c].n.begin(),alle[c].n.end(),vi[i]))
        vve.push_back(c);
  sort(vve.begin(),vve.end());
  vve.assign(vve.begin(),unique(vve.begin(),vve.end()));

  // create list of faces on the boundary from elements touching it
  // with appropriate number of nodes, flipping if necessary
  melem elem;
  ve.clear();
  for (unsigned i=0; i<vve.size(); ++i) {
    const vector< unsigned >& en = alle[vve[i]].n;
    elem.n.clear();
    vector< bool > nfound(en.size(),false);
    for (unsigned j=0; j<en.size(); ++j)
      if (nfound[j]=(count(vi.begin(),vi.end(),en[j])>0))
        elem.n.push_back(en[j]);

    if (elem.n.size()==2) {
      // line segment
      vector< unsigned >& bn = elem.n;
      const double v1x = vx[bn[1]]-vx[bn[0]],
                   v1y = vy[bn[1]]-vy[bn[0]],
                   v2x = (vx[en[0]]+vx[en[1]]+vx[en[2]])/3. - vx[bn[0]],
                   v2y = (vy[en[0]]+vy[en[1]]+vy[en[2]])/3. - vy[bn[0]];
      if (v1x*v2y-v1y*v2x<0.)
        swap(bn[0],bn[1]);
      ve.push_back(elem);
    }
    else if (elem.n.size()==3) {
      // triangle
      cerr << "error: triangle faces not implemented" << endl;
      throw 42;
    }

  }
}


// tecplot specific functions (input) ----------------------------------------


vector< string > f_plt::getVariables(const string& s)
{
  vector< string > r(1);

  // take out "VARIABLES" (up to '='), then take out '='
  const size_t eq(s.find('=',9));
  const string y = s.substr(eq!=std::string::npos? eq:11);
  istringstream ss(y);
  ss >> r.back();

  while (ss >> r.back()) {
    r.back() = trim(r.back(),"\"");
    r.push_back("");
  }
  r.pop_back();
  return r;
}


TecZone f_plt::getZoneHeader(const string& s, string& n, mtype& t)
{
  TecZone tz;

  istringstream ss;
  string key, value, dummy;  // dummy holds "="

  tz.isblock  = false;
  tz.isshared = false;

  vector< string > vs = splitstring(s);
  for (unsigned k=0; k<vs.size(); ++k) {
    ss.clear();  ss.str(vs[k]);  ss >> key;  key = upper(key);

    if      (key=="T")                                 tz.t = vs[k].substr(4);
    else if (key=="I" || key=="N" || key=="NODES")     ss >> dummy >> tz.i;
    else if (key=="J" || key=="E" || key=="ELEMENTS")  ss >> dummy >> tz.j;
    else if (key=="ZONETYPE" || key=="ET") {
      ss >> dummy >> value;
      value = upper(value);
      if      (value=="FELINESEG"       || value=="LINESEG")        { tz.nenodes = 2;  tz.et=FELINESEG;       }
      else if (value=="FETRIANGLE"      || value=="TRIANGLE")       { tz.nenodes = 3;  tz.et=FETRIANGLE;      }
      else if (value=="FEQUADRILATERAL" || value=="QUADRILATERAL")  { tz.nenodes = 4;  tz.et=FEQUADRILATERAL; }
      else if (value=="FETETRAHEDRON"   || value=="TETRAHEDRON")    { tz.nenodes = 4;  tz.et=FETETRAHEDRON;   }
      else if (value=="FEBRICK"         || value=="BRICK")          { tz.nenodes = 8;  tz.et=FEBRICK;         }
    }
    else if (key=="DATAPACKING" || key=="F") {
      ss >> dummy >> value;
      value = upper(value);
      tz.isblock = (value=="BLOCK");
    }
    else if (key=="VARSHARELIST") {
      tz.isshared = true;
    }

  }
  n = tz.t;
  t = tz.et;
  return tz;
}


void f_plt::readZoneNodeValues(ifstream& f, vector< vector< double > >& vv, unsigned N, unsigned Nvars, bool isblock)
{
  vv.resize(Nvars,vector< double >(N));
  if (isblock) {
    for (unsigned i=0; i<Nvars; ++i)
      for (unsigned l=0; l<N; ++l)
        f >> vv[i][l];
  }
  else {
    for (unsigned l=0; l<N; ++l)
      for (unsigned i=0; i<Nvars; ++i)
        f >> vv[i][l];
  }
  string s;  // make sure line is terminated
  getline(f,s);
}


void f_plt::readZoneConnectivity(ifstream& f, vector< melem >& ve, unsigned N, unsigned Nenodes, mtype& t)
{
  string s;
  unsigned n;
  melem e;  // default element
  e.n.assign(Nenodes,0);
  ve.assign(N,e);  // resize connectivity
  for (unsigned i=0; i<N; ++i) {

    // read indices
    vector< unsigned >& en = ve[i].n;
    for (unsigned j=0; j<Nenodes; ++j) {
      f >> n;
      en[j] = n-1;  // 0-based index
    }

    // correct numberings
    switch (t) {
      case PRISM3: {        // corrected from FEBRICK (see below)
        en.erase(en.begin()+3);
        en.resize(6);
      } break;
      case PYRAMID4: {      // corrected from FEBRICK (see below)
        swap(en[2],en[3]);
        en.resize(5);
      } break;
      case FEBRICK: {       // represents also PRISM3 and PYRAMID4
        if (i==0 && en[2]==en[3] && en[6]==en[7]) {
          t = PRISM3;
          en.erase(en.begin()+3);
          en.resize(6);
        }
        else if (i==0 && ((en[4]==en[5])==en[6])==en[7]) {
          t = PYRAMID4;
          en.resize(5);
        }
        else {
          swap(en[2],en[3]); swap(en[6],en[7]);
        }
      } break;
      case FEPOLYGON:       // don't know what to do
      case FEPOLYHEDRON:    // ...
      case ORDERED:         // nothing to do
      case FEQUADRILATERAL: // ...
      case FETETRAHEDRON:   // ...
      case FELINESEG:       // ...
      case FETRIANGLE:      // ...
      default: break;
    }
  }
  getline(f,s);  // make sure line is terminated
}


// tecplot specific functions (output) ---------------------------------------


string f_plt::setVariables(const vector< string >& vn)
{
  ostringstream ss;
  ss << "VARIABLES =";
  for (unsigned v=0; v<vn.size(); ++v)
    ss << " \"" << vn[v] << '\"';
  return ss.str();
}


string f_plt::setZoneHeader(const TecZone& tz, const unsigned Nvars, const double solutiontime)
{
  const string zt(tz.et==ORDERED?         "ORDERED"         :
                 (tz.et==FELINESEG?       "FELINESEG"       :
                 (tz.et==FETRIANGLE?      "FETRIANGLE"      :
                 (tz.et==FEQUADRILATERAL? "FEQUADRILATERAL" :
                 (tz.et==FETETRAHEDRON?   "FETETRAHEDRON"   :
                 (tz.et==FEBRICK?         "FEBRICK"         :
                 (tz.et==FEPOLYGON?       "FEPOLYGON"       :
                 (tz.et==FEPOLYHEDRON?    "FEPOLYHEDRON"    :
                                          "" ))))))));
  if (!zt.length())
    return string();

  ostringstream ss;
  ss << "ZONE"
     << " T=\"" << tz.t << "\""
     << " SOLUTIONTIME=" << solutiontime
     << " DATAPACKING=" << (tz.isblock? "BLOCK":"POINT");
  if (tz.isshared)
    ss << " VARSHARELIST=([1-" << Nvars << "]=1)";
  if (tz.j)
    ss << " N=" << tz.i << " E=" << tz.j << " ZONETYPE=" << zt;
  else
    ss << " I=" << tz.i;
  return ss.str();
}


void f_plt::writeZoneNodeValues(ofstream& f, const vector< vector< double > >& vv, bool isblock, const int& sharezone)
{
  if (sharezone>=0) {
    // zone picks values from somewhere else, do nothing
  }
  else if (isblock) {

    for (unsigned v=0; v<vv.size(); ++v)
      writeZoneNodeValues(f,vv[v]);

  }
  else {

    const unsigned Nvars = vv.size();
    const unsigned Nnode = (Nvars? vv[0].size():0);
    for (unsigned n=0; n<Nnode; ++n) { for (unsigned v=0; v<Nvars; ++v) {
        f << ' ' << vv[v][n];
    } f << endl; }

  }
}


void f_plt::writeZoneNodeValues(ofstream& f, const vector< double >& v)
{
  for (unsigned i=0; i<v.size(); ++i) {
    f << ' ' << v[i];
    if (!((i+1)%100))
      f << endl;  // long lines makes tecplot cough
  }
  f << endl;
}


void f_plt::writeZoneConnectivity(ofstream& f, const vector< melem >& ve, const mtype& t)
{
  switch (t) {
    case FELINESEG:
    case FETRIANGLE:
    case FEQUADRILATERAL:
    case FETETRAHEDRON:
    case FEPOLYGON:
    case FEPOLYHEDRON: {
      for (unsigned c=0; c<ve.size(); ++c) { for (unsigned i=0; i<ve[c].n.size(); ++i) {
          f << ' ' << ve[c].n[i]+1;  // 1-based index
      } f << endl; }
    } break;
    case FEBRICK: {
      for (unsigned c=0; c<ve.size(); ++c) {
        const vector< unsigned >& en = ve[c].n;
        f << ' ' << en[0]+1 << ' ' << en[1]+1 << ' ' << en[3]+1 << ' ' << en[2]+1
          << ' ' << en[4]+1 << ' ' << en[5]+1 << ' ' << en[7]+1 << ' ' << en[6]+1
          <<  endl;
      }
    } break;
    case ORDERED:    // no connectivity
    case PRISM3:     // not supported
    case PYRAMID4:   // not supported
    default: break;
  }
}


// string manipulation -------------------------------------------------------


string f_plt::trimright(const string& s, const string& t)
{
  string str = s;
  return str.erase(s.find_last_not_of(t)+1);
}


string f_plt::trimleft(const string& s, const string& t)
{
  string str = s;
  return str.erase(0,s.find_first_not_of(t));
}


string f_plt::trim(const string& s, const string& t)
{
  string str = s;
  return trimleft(trimright(str,t),t);
}


vector< string > f_plt::splitstring(const string& s)
{
  string y;  // a string easier to parse
  bool instr = false;  // if we are in a string ("...")
  bool inset = false;  // if we are in a set ((...))
  for (unsigned i=0; i<s.size(); ++i) {
    const char c = s[i];
    if      (c=='"')  instr=!instr;
    else if (c=='(')  inset=true;
    else if (c==')')  inset=false;

    if      (instr) { if (c!='"')  y+= c; }
    else if (inset) { if (c!=' ')  y+= c; }
    else if (c=='=')  y += " = ";
    else if (c==',')  y += ' ';
    else if (c!='"')  y += c;
  }

  vector< string > lines(1);
  string s1,
         s2;
  istringstream ss(y);
  while (ss >> s2) {
    if (s2=="=") {
      string& line = lines.back();
      line = trim(line," ");
      if (line.length())
        lines.push_back("");
    }
    lines.back() += ' ' + s1;
    s1 = s2;
  }
  lines.back() += ' ' + s1;
  string& line = lines.back();
  line = trim(line," ");

  return lines;
}


string f_plt::upper(const string& s)
{
  string r(s);
  std::transform(r.begin(),r.end(),r.begin(), (int(*)(int)) toupper);
  return r;
}

