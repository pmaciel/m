
#include <sstream>
#include "boost/progress.hpp"
#include "mfactory.h"
#include "t_map.h"

using namespace std;
using namespace m;


Register< mtransform,t_map > mt_map(4,"-tmap","[str] [str] [real] map 3D surface mesh onto another, given",
                                      "",     "  [str], (old) zone to map against, and",
                                      "",     "  [str], (new) surface mesh file name, and",
                                      "",     "  [real], (new) surface mesh scaling factor" );


void t_map::transform(GetPot& o, mmesh& m)
{
  if (m.d()!=3)
    return;
  const string o_mold_str = o.get(o.inc_cursor(),"");
  const string o_mnew_str = o.get(o.inc_cursor(),"");
  const double o_mnew_f   = o.get(o.inc_cursor(),0.);


  cout << "info: old mesh: setup (" << o_mold_str << ")..." << endl;
  // subset current mesh with just the intended zone
  mmesh mold = m.extract(o_mold_str);
  if (mold.d()!=3 || mold.d(0)!=2 || mold.z()!=1) {
    cerr << "error: mesh must satisfy: d=" << mold.d() << "==3, d(0)=" << mold.d(0) << "==2 and z()==" << mold.z() << "==1" << endl;
    throw 42;
  }
  cout << "info: old mesh: setup." << endl;


  cout << "info: old mesh: projection..." << endl;
  const mpoint          n = findnormal(mold.vz.front(),mold.vv);  // normal
  pair< mpoint,mpoint > p = findpoints(mold.vz.front(),mold.vv);  // center/far nodes
  const mpoint w = ((p.second-p.first)^n).normalize();
  const mpoint u = n^w;  // parametric base
  const mpoint v = n^u;  // ...
  for (unsigned i=0; i<mold.n(); ++i) {
    const mpoint a = mpoint(mold.vv[0][i],mold.vv[1][i],mold.vv[2][i]) - p.first;
    mold.vv[0][i] = a&u;
    mold.vv[1][i] = a&v;
    mold.vv[2][i] = 0.;
  }
  mold.vn.resize(3); mold.vn[0]="u"; mold.vn[1]="v"; mold.vn[2]="w";
  mold.vv.resize(3);
  const double mold_rmax = mpoint((p.second-p.first)&u,(p.second-p.first)&v,0.).length();
  cout << "info: old mesh: radius = " << mold_rmax << endl;
  cout << "info: old mesh: projection." << endl;


  cout << "info: old mesh: edge projection..." << endl;
  mmesh edge = m.edge(o_mold_str);
  for (unsigned i=0; i<edge.n(); ++i) {
    const mpoint a = mpoint(edge.vv[0][i],edge.vv[1][i],edge.vv[2][i]) - p.first;
    edge.vv[0][i] = a&u;
    edge.vv[1][i] = a&v;
    edge.vv[2][i] = 0.;
  }
  edge.vn.resize(2); edge.vn[0]="u"; edge.vn[1]="v";
  edge.vv.resize(2);
  cout << "info: old mesh: edge projection." << endl;


  cout << "info: new mesh: setup (" << o_mnew_str << ")..." << endl;
  mmesh mnew;
  {
    // get a mesh reader and make it work
    // (including a stupid way to set GetPot, should only be an input filename)
    int argc = 2;
    char *argv[2];
    argv[0] = new char[1];              argv[0][0] = '\0';
    argv[1] = new char[o_mnew_str.length()+1];  strcpy(argv[1],o_mnew_str.c_str());
    GetPot o(argc,argv);
    const string::size_type idx = o_mnew_str.rfind('.');
    const string ext(o_mnew_str.substr(idx!=string::npos? idx:0));
    mfactory< mfinput >* fi = mfactory< mfinput >::instance();
    if (fi->search(ext.c_str())) {
      mfinput* p = fi->Create(ext);
      p->read(o,mnew);
      delete p;
    }
    delete[] argv[0];
    delete[] argv[1];

    // set size factor
    for (unsigned i=0; i<mnew.v(); ++i)
      for (vector< double >::iterator v=mnew.vv[i].begin(); v!=mnew.vv[i].end(); ++v)
        *v *= o_mnew_f;

    if (mnew.d()==2 && mnew.d(0)==2) {
      mnew.vn.resize(3);  mnew.vn[0]="u";  mnew.vn[1]="v";  mnew.vn[2]="w";
      mnew.vv.resize(3);  mnew.vv[2].assign(mnew.n(),0.);
    }
    if (mnew.d()!=3 || mnew.d(0)!=2 || mnew.z()!=1) {
      cerr << "error: mesh must satisfy: d=" << mnew.d() << "==3, d(0)=" << mnew.d(0) << "==2 and z()==" << mnew.z() << "==1" << endl;
      throw 42;
    }

    // clipping
    const double mnew_rclip = mold_rmax + 2*o_mnew_f;
    vector< melem > e2n;
    e2n.swap(mnew.vz[0].e2n);
    const vector< double >& u = mnew.vv[0];
    const vector< double >& v = mnew.vv[1];
    for (vector< melem >::const_iterator e=e2n.begin();e!=e2n.end(); ++e)
      for (vector< unsigned >::const_iterator n=(*e).n.begin();n!=(*e).n.end(); ++n)
        if (mpoint(u[*n],v[*n],0.).length()<mnew_rclip) {
          mnew.vz[0].e2n.push_back(*e);
          break;
        }
  }
  mnew.compress();
  cout << "info: new mesh: setup." << endl;


  cout << "info: new mesh: division by sampling..." << endl;
  {
    const unsigned Ntsample = 100;//200;   // samples per triangle
    const unsigned Nlsample = 1000;//2000;  // samples per line segment
    mpoint pts[3];
    vector< mpoint > sspt;  sspt.reserve(Ntsample*mold.e(0));  // points in triangles
    vector< mpoint > sspl;  sspl.reserve(Nlsample*edge.e(0));  // points in line segments
    for (vector< melem >::const_iterator e=mold.vz[0].e2n.begin(); e!=mold.vz[0].e2n.end(); ++e) {
      pts[0] = mpoint( mold.vv[0][e->n[0]],mold.vv[1][e->n[0]],0. );
      pts[1] = mpoint( mold.vv[0][e->n[1]],mold.vv[1][e->n[1]],0. );
      pts[2] = mpoint( mold.vv[0][e->n[2]],mold.vv[1][e->n[2]],0. );
      for (unsigned j=0; j<Ntsample; ++j)
        sspt.push_back(intriangle(pts[0],pts[1],pts[2]));
    }
    for (vector< melem >::const_iterator e=edge.vz[0].e2n.begin(); e!=edge.vz[0].e2n.end(); ++e) {
      pts[0] = mpoint( edge.vv[0][e->n[0]],edge.vv[1][e->n[0]],0. );
      pts[1] = mpoint( edge.vv[0][e->n[1]],edge.vv[1][e->n[1]],0. );
      for (unsigned j=0; j<Nlsample; ++j)
        sspl.push_back(inlinesegment(pts[0],pts[1]));
    }

    vector< melem > e2n_t;
    vector< melem > e2n_l;
    boost::progress_display pbar(mnew.vz[0].e2n.size());
    for (vector< melem >::const_iterator e=mnew.vz[0].e2n.begin(); e!=mnew.vz[0].e2n.end(); ++e, ++pbar) {
      pts[0] = mpoint( mnew.vv[0][e->n[0]],mnew.vv[1][e->n[0]],0. );
      pts[1] = mpoint( mnew.vv[0][e->n[1]],mnew.vv[1][e->n[1]],0. );
      pts[2] = mpoint( mnew.vv[0][e->n[2]],mnew.vv[1][e->n[2]],0. );
      bool in = false;
      for (vector< mpoint >::const_iterator p=sspl.begin(); p!=sspl.end() && !in; ++p)
        if (p->intriangle(pts[0],pts[1],pts[2])) {
          e2n_l.push_back(*e);
          in = true;
        }
      for (vector< mpoint >::const_iterator p=sspt.begin(); p!=sspt.end() && !in; ++p)
        if (p->intriangle(pts[0],pts[1],pts[2])) {
          e2n_t.push_back(*e);
          in = true;
        }
    }

    mnew.vz.resize(2);
    e2n_t.swap(mnew.vz[0].e2n);  // this part is alright
    e2n_l.swap(mnew.vz[1].e2n);  // this part needs remeshing
  }
  mnew.compress();
  cout << "info: new mesh: division by sampling." << endl;


  cout << "info: mesh edge adaptation..." << endl;
  {
    const unsigned Ntri = mnew.e(1);
    const unsigned Npnt = edge.n();

    vector< vector< mpoint > > vtriangles(Ntri,vector< mpoint >(3));
    for (unsigned i=0; i<Ntri; ++i) {
      const vector< unsigned >& en = mnew.vz[1].e2n[i].n;
      vtriangles[i][0] = mpoint( mnew.vv[0][en[0]],mnew.vv[1][en[0]],0. );
      vtriangles[i][1] = mpoint( mnew.vv[0][en[1]],mnew.vv[1][en[1]],0. );
      vtriangles[i][2] = mpoint( mnew.vv[0][en[2]],mnew.vv[1][en[2]],0. );
    }
    vector< mpoint > vpoints(Npnt);
    for (unsigned i=0; i<Npnt; ++i)
      vpoints[i] = mpoint( edge.vv[0][i],edge.vv[1][i],0. );

    vector< vector < bool > > contains(Ntri,vector< bool >(Npnt,false));
    {
      unsigned it, ip;
      vector<      mpoint      >::const_iterator p;
      vector< vector< mpoint > >::const_iterator t;
      for (p=vpoints.begin(), ip=0; p!=vpoints.end(); ++p, ++ip)
        for (t=vtriangles.begin(), it=0; t!=vtriangles.end(); ++t, ++it)
          contains[it][ip] = p->intriangle((*t)[0],(*t)[1],(*t)[2]);
    }

    //vector< pair< mpoint, vector< mpoint > > > adapt;
    vector< unsigned > c(3,0);
    for (vector< vector< bool > >::const_iterator t=contains.begin(); t!=contains.end(); ++t)
      ++c[ min((unsigned) 2,(unsigned) count(t->begin(),t->end(),true)) ];
    cout << "info: 0/3 adaptation triangles with 0/1/>1 edge points: " << c[0] << '/' << c[1] << '/' << c[2] << endl;

    c.assign(3,0);
    for (vector< vector< bool > >::const_iterator t=contains.begin(); t!=contains.end(); ++t)
      ++c[ min((unsigned) 2,(unsigned) count(t->begin(),t->end(),true)) ];
    cout << "info: 1/3 adaptation triangles with 0/1/>1 edge points: " << c[0] << '/' << c[1] << '/' << c[2] << endl;

  }
  cout << "info: mesh edge adaptation." << endl;







#if 1
static unsigned notthefirst = 0;
ofstream f3d("dump.3d.plt",notthefirst? ios_base::app:ios_base::trunc);
ofstream f2d("dump.2d.plt",notthefirst? ios_base::app:ios_base::trunc);
if (!(notthefirst++)) {
f3d << "VARIABLES = x y z" << endl;
f2d << "VARIABLES = u v" << endl;
}
dump2d(f2d,edge.vz.front(),edge.vv);
dump2d(f2d,mnew.vz.back(),mnew.vv);
#endif

}


pair< mpoint,mpoint > t_map::findpoints(const mzone& z, const vector< vector< double > >& vv)
{
  pair< mpoint,mpoint > r;
  mpoint& pc = r.first;
  mpoint& pf = r.second;

  // find center node as average of all coordinates
  const unsigned Nnode = vv[0].size();
  for (unsigned i=0; i<Nnode; ++i)
    pc += mpoint(vv[0][i],vv[1][i],vv[2][i]);
  pc /= (double) Nnode;

  // find farthest from center node index
  double   max_d = 0.;  // max. distance (squared)
  unsigned max_i = 0;   // ...  node index
  for (unsigned i=0; i<vv[0].size(); ++i) {
    const double d = pc.norm2(mpoint(vv[0][i],vv[1][i],vv[2][i]));
    if (d>max_d) {
      max_d = d;
      max_i = i;
    }
  }
  pf = mpoint(vv[0][max_i],vv[1][max_i],vv[2][max_i]);

  return r;
}


mpoint t_map::findnormal(const mzone& z, const vector< vector< double > >& vv)
{
  // get (and distribute) element normal
  mpoint cn;
  for (unsigned i=0; i<z.e2n.size(); ++i) {
    const vector< unsigned >& en = z.e2n[i].n;
    const mpoint p0 = mpoint(vv[0][en[0]],vv[1][en[0]],vv[2][en[0]]);
    const mpoint p1 = mpoint(vv[0][en[1]],vv[1][en[1]],vv[2][en[1]]);
    const mpoint p2 = mpoint(vv[0][en[2]],vv[1][en[2]],vv[2][en[2]]);
    cn += (p2-p0)^(p1-p0);
  }
   return cn.normalize();
}


void t_map::dump2d(ostream& out, const mzone& z, const vector< vector< double > >& vv)
{
  const unsigned Nnode = (vv.size()? vv[0].size():0);
  const unsigned Nelem  = z.e2n.size();
  const unsigned Nnelem = (Nelem? z.e2n[0].n.size():0);
  out << "ZONE T=\"" << z.n << "\" "
      << " ZONETYPE=" << (Nnelem==3? "FETRIANGLE":
                         (Nnelem==2? "FELINESEG":
                                     "ORDERED" ))
      << " DATAPACKING=POINT N=" << Nnode << " E=" << Nelem << endl;
  for (unsigned i=0; i<Nnode; ++i)
    out << vv[0][i] << ' ' << vv[1][i] << endl;
  for (unsigned i=0; i<Nelem; ++i) {
    for (vector< unsigned >::const_iterator n=z.e2n[i].n.begin(); n!=z.e2n[i].n.end(); ++n)
      out << ' ' << *n+1;
    out << endl;
  }
}


void t_map::dump3d(ostream& out, const mzone& z, const vector< vector< double > >& vv)
{
  const unsigned Nnode = (vv.size()? vv[0].size():0);
  const unsigned Nelem  = z.e2n.size();
  const unsigned Nnelem = (Nelem? z.e2n[0].n.size():0);
  out << "ZONE T=\"" << z.n << "\" "
      << " ZONETYPE=" << (Nnelem==3? "FETRIANGLE":
                         (Nnelem==2? "FELINESEG":
                                     "ORDERED" ))
      << " DATAPACKING=POINT N=" << Nnode << " E=" << Nelem << endl;
  for (unsigned i=0; i<Nnode; ++i)
    out << vv[0][i] << ' ' << vv[1][i] << ' ' << vv[2][i] << endl;
  for (unsigned i=0; i<Nelem; ++i) {
    for (vector< unsigned >::const_iterator n=z.e2n[i].n.begin(); n!=z.e2n[i].n.end(); ++n)
      out << ' ' << *n+1;
    out << endl;
  }
}


void t_map::dumpuv(ostream& out, const mzone& z, const vector< vector< double > >& vv, const mpoint& c, const mpoint& u, const mpoint& v)
{
  const unsigned Nnode = (vv.size()? vv[0].size():0);
  const unsigned Nelem  = z.e2n.size();
  const unsigned Nnelem = (Nelem? z.e2n[0].n.size():0);
  out << "ZONE T=\"" << z.n << "\" "
      << " ZONETYPE=" << (Nnelem==3? "FETRIANGLE":
                         (Nnelem==2? "FELINESEG":
                                     "ORDERED" ))
      << " DATAPACKING=POINT N=" << Nnode << " E=" << Nelem << endl;
  for (unsigned i=0; i<Nnode; ++i) {
    const mpoint a(mpoint(vv[0][i],vv[1][i],vv[2][i]) - c);
    out << (a&u) << ' ' << (a&v) << endl;
  }
  for (unsigned i=0; i<Nelem; ++i) {
    for (vector< unsigned >::const_iterator n=z.e2n[i].n.begin(); n!=z.e2n[i].n.end(); ++n)
      out << ' ' << *n+1;
    out << endl;
  }
}

