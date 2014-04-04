
#include <sstream>
#include "mfactory.h"
#include "f_miotras.h"

using namespace std;
using namespace m;


Register< mfinput,f_miotras >  mf_miotras1(2,".grid","Miotras input format (mesh)",
                                             ".flow","Miotras input format (mesh and solution)");
Register< mfoutput,f_miotras > mf_miotras2(4,".miotras","Miotras output format",
                                             "","--miotras-boundaries [i:i:...] output boundary types",
                                             "","  0:insulator, 1:electrode, 2:wall, 3:inlet, 4:outlet",
                                             "","(must specify for all boundaries, default: 0:0:...)");


// useful definitions for this format
enum BCType {
  NONE = 0, INSULATOR = 1, ELECTRODE = 2,
  WALL = 3, INLET     = 4, OUTLET    = 5 };

const unsigned BCTypePriority[] = {
   0,  10,  20,
  30,  40,  50 };


void f_miotras::read(GetPot& o, mmesh& m)
{
  // try to access file
  const string fn(o.get(o.inc_cursor(),""));
  ifstream f(fn.c_str());
  if (!f) {
    cerr << "error accessing file: \"" << fn << "\"." << endl;
    throw 42;
  }

  // detect line ending type, read first line to detect file type, and choose
  // appropriate reading method
  string
    line,
    line_cr;
  getline(f,line,   '\n');  f.seekg(0,ios::beg);
  getline(f,line_cr,'\r');  f.seekg(0,ios::beg);
  if (line.length()>line_cr.length())
    line = line_cr;

  if      (line=="Miotras Genereated Grid") read_grid(f,m);
  else if (line=="Fluid flow witk k-omega") read_flow(f,m);
  else {
    cerr << "error: the provided file: \"" << fn << "\" is not of the expected format." << endl;
    throw 42;
  }

  // close file
  f.close();
}


void f_miotras::write(GetPot& o, const mmesh& m)
{
  const string fn(o.get(o.inc_cursor(),""));
  const unsigned dim = m.d();
  if (dim!=2) {
    cerr << "number of dimensions should be 2!" << endl;
    throw 42;
  }
  string btype(o.follow("","--miotras-boundaries"));


  // -- reorder zones, setting the first one as the volume

  vector< const mzone* > z;
  for (unsigned i=0; i<m.vz.size(); ++i) {
    if (m.d(i)==dim && m.vz[i].t==FETRIANGLE) {
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

  const unsigned Nnode = (unsigned) m.vv[0].size();
  const unsigned Nbc = (unsigned) z.size()-1;


  // -- open files --

  ofstream fgrid((fn+".grid").c_str(),ios_base::trunc);
  ofstream fflow((fn+".flow").c_str(),ios_base::trunc);
  ofstream fdat((fn+".dat").c_str(),ios_base::trunc);
  if (!fgrid) {
    cerr << "error accessing file \"" << (fn+".grid") << "\"" << endl;
    throw 42;
  }
  if (!fflow) {
    cerr << "error accessing file \"" << (fn+".flow") << "\"" << endl;
    throw 42;
  }
  if (!fdat) {
    cerr << "error accessing file \"" << (fn+".dat") << "\"!" << endl;
    throw 42;
  }
  fgrid.precision(15);
  fflow.precision(15);
  fdat.precision(15);


  // -- .grid: header --

  fgrid << "Miotras Genereated Grid" << endl
        << "2D" << endl
        << "6.1" << endl
        << "1" << endl
        << Nnode << endl;
  fflow << "Fluid flow witk k-omega" << endl
        << Nnode << endl;


  // -- .grid: nodal information --

  // recognize boundary types
  vector< BCType > ztypes(1+Nbc,NONE);
  if (Nbc==(unsigned) count(btype.begin(),btype.end(),',')+1) {
    replace(btype.begin(),btype.end(),',',' ');
    istringstream is(btype);
    for (unsigned t=1; t<ztypes.size(); ++t) {
      int T;
      is >> T;
      ztypes[t] = (BCType) (T+1);
    }
  }

  // output boundary types
  for (unsigned t=1; t<=Nbc; ++t) {
    string n = "<unknown>";
    switch (ztypes[t]) {
      case NONE:       n = "none";       break;
      case INSULATOR:  n = "insulator";  break;
      case ELECTRODE:  n = "electrode";  break;
      case WALL:       n = "wall";       break;
      case INLET:      n = "inlet";      break;
      case OUTLET:     n = "outlet";     break;
    }
    cout << "boundary " << t << "/" << Nbc << ": " << n << endl;
  }

  // boundary nodes belong to and nodal coordinates/values
  {
    vector< BCType >   vnt(Nnode,NONE);  // boundary types where nodes sit
    vector< unsigned > vni(Nnode,0);     // ... indices where nodes sit

    for (unsigned t=1; t<=Nbc; ++t) {
      for (unsigned c=0; c<z[t]->e2n.size(); ++c)
        for (unsigned i=0; i<z[t]->e2n[c].n.size(); ++i) {
          const unsigned N = z[t]->e2n[c].n[i];
          if ( BCTypePriority[ vnt[N] ] < BCTypePriority[ ztypes[t] ]) {
            vnt[N] = ztypes[t];
            vni[N] = t;
          }
        }
    }

    for (unsigned n=0; n<Nnode; ++n) {
      fgrid << vni[n] << endl;
      fflow << vni[n] << "\t"
            << m.vv[0][n] << "\t" << m.vv[1][n] << "\t"  // x,y
            << "0\t0\t0" << "\t"                         // p,u,v
            << "0\t0" << endl;                           // k,w
    }

    // coordinates
    for (unsigned n=0; n<Nnode; ++n)
      fgrid << m.vv[0][n] << "\t" << m.vv[1][n] << endl;
  }


  // -- .grid: elements information --

  unsigned Nbe = 0;
  for (unsigned t=1; t<=Nbc; ++t)
    Nbe += z[t]->e2n.size();
  fgrid << z[0]->e2n.size()+Nbe << "\t" << z[0]->e2n.size() << "\t" << Nbe << endl;
  fflow << z[0]->e2n.size()+Nbe << " "  << z[0]->e2n.size() << " "  << Nbe << endl;

  for (unsigned c=0; c<z[0]->e2n.size(); ++c)
    fgrid << 0 << endl;
  for (unsigned t=1; t<=Nbc; ++t) {
    for (unsigned c=0; c<z[t]->e2n.size(); ++c)
      fgrid << t << endl;
  }

  {
    // duplicate connectivity to include boundaries
    // (to search for boundary element indices)
    vector< melem > ve(z[0]->e2n);
    ve.reserve(ve.size()+Nbe);
    for (unsigned t=1; t<=Nbc; ++t)
      for (unsigned c=0; c<z[t]->e2n.size(); ++c)
        ve.push_back(z[t]->e2n[c]);

    // mark elements according to boundary indices
    vector< int > vbzones(ve.size());
    {
      vector< int >::iterator i = vbzones.begin();
      for (unsigned c=0; c<z[0]->e2n.size(); ++c, ++i)
        *i = -1;
      for (unsigned t=1; t<=Nbc; ++t)
        for (unsigned c=0; c<z[t]->e2n.size(); ++c, ++i)
          *i = (int) t;
    }

    // output volume and boundaries connectivity (1-based index)
    for (unsigned c=0; c<ve.size(); ++c) {
      if (!(c%1000))
        cout << "finding elements neighbours: "
             << (int)(100.*((double)c)/((double)ve.size())) << "%" << endl;

      fgrid << 1          << "\t" << ve[c].n.size() << "\t";
      fflow << vbzones[c] << " "  << ve[c].n.size() << "\t";

      if (ve[c].n.size()==2) {
        // for some reason boundary elements have to be reversed
        fgrid << ve[c].n[1]+1 << "\t" << ve[c].n[0]+1 << "\t";
        fflow << ve[c].n[1]+1 << " "  << ve[c].n[0]+1 << " ";
      }
      else {
        for (unsigned i=0; i<ve[c].n.size(); ++i) {
          fgrid << ve[c].n[i]+1 << "\t";
          fflow << ve[c].n[i]+1 << " ";
        }
      }

      vector< unsigned > vneigh = findneighbours(ve,ve[c].n);
      fgrid <<         vneigh.size() << "\t";
      fflow << "\t" << vneigh.size() << " ";
      for (unsigned i=0; i<vneigh.size(); ++i) {
        fgrid << "1\t1\t" << vneigh[i]+1 << "\t";
        fflow << "1 1 "   << vneigh[i]+1 << " ";
      }

      fgrid << "0 0 0" << endl;
      fflow <<            endl;

    }
    if ((ve.size()%1000))
      cout << "finding elements neighbours: 100%" << endl;
  }


  // -- .dat --
  fdat << "Miotrasfile version 1.0" << endl  // FileString
       << "1" << endl                 // numbDomains
       << "1" << endl                 // thirdDimension
       << "1" << endl                 // readDimension probDim (1:TWOD;2:AX)
       << "4" << endl                 // readUnits unit (1:MM;2:CM;3:DM;4:M)
       << "0" << endl                 // numbparams
       << "0" << endl;                // size (none, so no size lines)
  fdat << Nbc << endl;
  for (unsigned t=1; t<=Nbc; ++t)
    fdat << "1\t0\t"                  // numblines lines with domain; type
         << "9\t9\t0\t"               // first; second and mid points
         << z[t]->e2n.size() << "\t"  // numbofelem
         << "0\t0\t"                  // distribmid; elemdistrib
         << ztypes[t] << "\t"         // bctype
         << "0\t1e-06\t"              // layerpresent; firstlayer
         << "<potimposed>\t<potential>\t"
         << "<currentimposed>\t<current>\t"
         << "0\t0\t0.001\t1\t"              // resistive; resistivity; thickness; contactbegin
         << "0\t0\t0\t0\t0\t0\t0" << endl;  // 3 times (val, imp); ElecType,if (Type == 2): CompShapeName; position; position;
  fdat << "<ionsystemindex>" << endl
       << "<nstotiterations>" << endl
       << "<dummy>" << endl
       << "<dummy>" << endl
       << "<elchemmaxiterations>" << endl
       << "<elchemcfl>" << endl
       << "<elchemconvergencefactor>" << endl;


  // -- close files --

  fgrid.close();
  fflow.close();
  fdat.close();
}


void f_miotras::read_flow(std::ifstream& f, mmesh& m)
{
  unsigned dummy, Nnode, Nelem, Nbelem;
  string s;

  // skip initial line
  getline(f,s);

  // allocate memory for node coordinates and solution fields
  {
    getline(f,s);
    istringstream ss(s);
    ss >> Nnode;
    m.vv.assign(7,vector< double >(Nnode,0.));
    m.vn.resize(7);
    m.vn[0] = "x";
    m.vn[1] = "y";
    m.vn[2] = "vx";
    m.vn[3] = "vy";
    m.vn[4] = "p";
    m.vn[5] = "k";
    m.vn[6] = "omega";
  }

  // read node coordinates and solution field
  for (unsigned i=0; i<Nnode; ++i) {
    getline(f,s);
    istringstream ss(s);
    ss >> dummy;
    for (unsigned d=0; d<7; ++d)
      ss >> m.vv[d][i];
  }

  // set volume/boundary elements zone map and connectivity lists
  vector< unsigned > vezonemap;
  vector< vector< unsigned > > veconn;
  {
    getline(f,s);
    istringstream ss1(s);
    ss1 >> dummy >> Nelem >> Nbelem;

    vezonemap.assign(Nelem+Nbelem,0);
    veconn.resize(Nelem+Nbelem);
    for (unsigned i=0; i<Nelem+Nbelem; ++i) {
      getline(f,s);
      istringstream ss2(s);
      ss2 >> vezonemap[i];

      unsigned Nenode, n;
      ss2 >> Nenode;
      if (Nenode==2) {
        veconn[i].resize(2);
        ss2 >> n; veconn[i][1] = n-1;
        ss2 >> n; veconn[i][0] = n-1;
      }
      else if (Nenode==3) {
        veconn[i].resize(3);
        ss2 >> n; veconn[i][0] = n-1;
        ss2 >> n; veconn[i][1] = n-1;
        ss2 >> n; veconn[i][2] = n-1;
      }
    }
  }

  // allocate memory for zones
  m.vz.resize(1+*max_element(vezonemap.begin(),vezonemap.end()));
  melem e;
  m.vz[0].n = "InnerCells";
  m.vz[0].e2n.reserve(Nelem);
  m.vz[0].t = FETRIANGLE;
  for (unsigned t=1; t<m.vz.size(); ++t) {
    ostringstream name;
    name << "boundary" << t;
    e.n.assign(2,0);
    m.vz[t].n = name.str();
    m.vz[t].e2n.reserve(count(vezonemap.begin(),vezonemap.end(),t));
    m.vz[t].t = FELINESEG;
  }

  // fill elements connectivity
  for (unsigned ie=0; ie<Nelem+Nbelem; ++ie) {
    m.vz[ vezonemap[ie] ].e2n.push_back(e);
    m.vz[ vezonemap[ie] ].e2n.back().n.swap( veconn[ie] );
  }
}


void f_miotras::read_grid(std::ifstream& f, mmesh& m)
{
  unsigned dummy, Nnode, Nelem, Nbelem;
  string s;

  // skip 4 initial lines
  for (int i=0; i<4; ++i)
    getline(f,s);

  // allocate memory for node coordinates
  {
    getline(f,s);
    istringstream ss(s);
    ss >> Nnode;
    m.vv.assign(2,vector< double >(Nnode,0.));
    m.vn.resize(2);
    m.vn[0] = "x";
    m.vn[1] = "y";
  }

  // skip node zone flags
  for (unsigned i=0; i<Nnode; ++i)
    getline(f,s);

  // read node coordinates
  for (unsigned i=0; i<Nnode; ++i) {
    getline(f,s);
    istringstream ss(s);
    for (unsigned d=0; d<2; ++d)
      ss >> m.vv[d][i];
  }

  // set volume/boundary elements zone map
  vector< unsigned > vezonemap;
  {
    getline(f,s);
    istringstream ss1(s);
    ss1 >> dummy >> Nelem >> Nbelem;

    vezonemap.assign(Nelem+Nbelem,0);
    for (unsigned i=0; i<Nelem+Nbelem; ++i) {
      getline(f,s);
      istringstream ss2(s);
      ss2 >> vezonemap[i];
    }
  }

  // read elements flags and allocate memory for zones
  m.vz.resize(1+*max_element(vezonemap.begin(),vezonemap.end()));
  melem e;
  e.n.assign(3,0);
  m.vz[0].n = "InnerCells";
  m.vz[0].e2n.assign(Nelem,e);
  m.vz[0].t = FETRIANGLE;
  for (unsigned t=1; t<m.vz.size(); ++t) {
    ostringstream name;
    name << "boundary" << t;
    e.n.assign(2,0);
    m.vz[t].n = name.str();
    m.vz[t].e2n.assign(count(vezonemap.begin(),vezonemap.end(),t),e);
    m.vz[t].t = FELINESEG;
  }

  // read elements connectivity
  vector< unsigned > vccounter(m.vz.size(),0); // element counter, per zone
  for (unsigned i=0; i<vezonemap.size(); ++i) {
    const unsigned c = vccounter[ vezonemap[i] ]++;
    vector< unsigned >& en = m.vz[ vezonemap[i] ].e2n[ c ].n;
    getline(f,s);
    istringstream ss(s);
    if (en.size()==3)
      ss >> dummy >> dummy >> en[0] >> en[1] >> en[2];
    else if (en.size()==2)
      ss >> dummy >> dummy >> en[1] >> en[0];
    for (unsigned j=0; j<en.size(); ++j)
      --en[j];  // 0-based index
  }
  vezonemap.clear();
}


vector< unsigned > f_miotras::findneighbours(const vector< melem >& ve, const vector< unsigned >& n)
{
  if (n.size()==2) {

    // special case for boundary elements: returns first element
    // that shares 2 nodes, but has more than 2 nodes itself (inner)
    for (unsigned c=0; c<ve.size(); ++c) {
      if ( count(ve[c].n.begin(),ve[c].n.end(),n[0]) &&
           count(ve[c].n.begin(),ve[c].n.end(),n[1]) &&
           ve[c].n.size()>2 )
        return vector< unsigned >(1,c);
    }

  }
  else if (n.size()==3) {

    // find elements sharing a face (1 shared node only touches, 3 is the same
    // element) and set neighbour position based on the edges 1-2, 2-3 and 3-1
    vector< unsigned > vi(3,0);
    unsigned nnfound = 0;
    for (unsigned c=0; c<ve.size(); ++c) {
      const vector< unsigned >& enodes = ve[c].n;
      if ( (count(enodes.begin(),enodes.end(),n[0])? 1:0) +
           (count(enodes.begin(),enodes.end(),n[1])? 1:0) +
           (count(enodes.begin(),enodes.end(),n[2])? 1:0) == 2) {

        for (unsigned e=0; e<3; ++e) {
          vector< unsigned > edgenodes(2,n[e]);
                             edgenodes[1] = n[(e+1)%3];
          if ( count(enodes.begin(),enodes.end(),edgenodes[0]) &&
               count(enodes.begin(),enodes.end(),edgenodes[1]) ) {
            vi[e] = c;
            break;
          }
        }

        if (++nnfound==3)
          return vi;

      }
    }
    return vi;  // this is just for safety

  }

  return vector< unsigned >();
}


