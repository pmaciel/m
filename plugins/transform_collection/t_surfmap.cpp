
#include <set>
#include <numeric>
#include "boost/progress.hpp"
#include "ext/Vec.h"
#include "mfactory.h"
#include "mpoint.h"
#include "t_surfmap.h"

using namespace std;
using namespace m;


Register< mtransform,t_surfmap > mt_surfmap(5,"-tsurfmap","[str] [real] [real] [real] map a structured grid on a 3D surface mesh, given",
                                              "",         "with [str] zone to map against, and",
                                              "",         "with [real] structured grid u vertex distance",
                                              "",         "with [real] structured grid v vertex distance",
                                              "",         "with [real] uv projection angle" );


// useful definitions
#define NEXT(i) ((i)<2 ? (i)+1 : (i)-2)
#define PREV(i) ((i)>0 ? (i)-1 : (i)+2)
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

// edge definition

// triangle index/weight pair, edge, point and triangle definitions
typedef pair< int,double > WPAIR;
struct WPAIRLT {
  bool operator()(const WPAIR& p1, const WPAIR& p2) const { return p1.first<p2.first; }
};
typedef pair< int,int > EDGE;
struct EDGELT {
  bool operator()(const EDGE& e1, const EDGE& e2) const {
    return (e1.first<e2.first?   true  :
            e1.first>e2.first?   false :
            e1.second<e2.second? true  :
                                 false );
  }
};
typedef Vec< 3,double > PNT;
struct TRI {
  Vec< 3,int > v;  // (splitted) triangle mmesh node indices
  double       w;  // (splitted) tri. contribution weight (start with 1.)
  int          i;  // (original) tri. mmesh cell index
  double       a;  // (original) tri. area
};
inline ostream& operator<< (ostream& os, const TRI& t) {
  for (int i=0; i<3; ++i)
    os << ' ' << t.v[i]+1;
  return os << ' ';
}
struct MNODE {
  PNT nucsite_xyz;      // nucleation site (in xyz)
  PNT nucsite_uv;       // nucleation site (in uv)
  PNT    qcenter;       // quad. center (in uv)
  double qarea;         // quad. area (intersected with zone)
  vector< int    > ni;  // contributing triangle index
  vector< double > nw;  // contributing triangle weight
};


void t_surfmap::transform(GetPot& o, mmesh& m)
{
  // options
  const string zonename = o.get(o.inc_cursor(),"");
  const double du       = o.get(o.inc_cursor(),1.);
  const double dv       = o.get(o.inc_cursor(),1.);
  const double angled   = o.get(o.inc_cursor(),0.);


  cout << "info: setup projection vectors (" << zonename << ")..." << endl;
  // subset current mesh with just the intended zone
  mmesh mold = m.extract(zonename);
  if (mold.d()!=3 || mold.d(0)!=2 || mold.z()!=1) {
    cerr << "error: mesh must satisfy: d=" << mold.d() << "==3, d(0)=" << mold.d(0) << "==2 and z()==" << mold.z() << "==1" << endl;
    throw 42;
  }

  // set points (3d/uvw) and triangles, together with pr. center and normal
  vector< pair< PNT,PNT > > points(mold.n());
  vector< TRI > tri(mold.e(0));
  PNT pu,            // proj. vector u
      pv,            // ...          v
      pc(0.,0.,0.),  // proj. surface center
      pn(0.,0.,0.);  // ...           normal
  {
    for (unsigned i=0; i<mold.n(); ++i) {
      points[i].first = PNT(mold.vv[0][i],mold.vv[1][i],mold.vv[2][i]);
      pc += points[i].first;
    }
    pc /= (double) mold.n();

    for (unsigned i=0; i<mold.e(0); ++i) {
      const vector< unsigned >& en = mold.vz[0].e2n[i].n;
      const PNT tnormal = trinorm(points[en[0]].first,points[en[1]].first,points[en[2]].first);
      tri[i].v = Vec<3,int>(en[0],en[1],en[2]);
      tri[i].w = 1.;
      tri[i].i = (int) i;       // kept constant throughout splitting
      tri[i].a = len(tnormal);  // ...
      pn += tnormal;
    }
    pn /= (double) mold.e(0);

    // set projection vectors
    pv = pn % (points.front().first-pc);
    pu = pv % pn;
    normalize(pu);
    normalize(pv);

    // rotate projection vectors (don't rotate, doesn't work...)
    /*
    const double r[4] = { cos(angled*M_PI/180.), -sin(angled*M_PI/180.),
                          sin(angled*M_PI/180.),  cos(angled*M_PI/180.) };
    PNT pu2(pu);  pu2[0]=r[0]*pu[0]+r[1]*pu[1];  pu2[1]=r[2]*pu[0]+r[3]*pu[1];  pu=pu2;
    PNT pv2(pv);  pv2[0]=r[0]*pv[0]+r[1]*pv[1];  pv2[1]=r[2]*pv[0]+r[3]*pv[1];  pv=pv2;
    */

    // information
    cout << "info: projection vector u: " << pu << endl;
    cout << "info: projection vector v: " << pv << endl;
    cout << "info: projection vector n: " << pn << endl;
  }
  cout << "info: setup projection vectors." << endl;


  // setup an edge and mark its nodes
  set< EDGE,EDGELT > edge;
  for (unsigned i=0; i<tri.size(); ++i) {
    for (int e=0; e<3; ++e) {
      const int v1 = tri[i].v[e];
      const int v2 = tri[i].v[NEXT(e)];
      // opposing edges cancel each other
      if (!edge.erase(EDGE(v2,v1)))
        edge.insert(EDGE(v1,v2));
    }
  }
  vector< bool > isedge(mold.n(),false);
  for (set< EDGE,EDGELT >::const_iterator e=edge.begin(); e!=edge.end(); ++e)
    isedge[ e->first  ] = isedge[ e->second ] = true;


  cout << "info: project at rotation angle [deg]: " << angled << "..." << endl;
  const double r[4] = { cos(angled*M_PI/180.), -sin(angled*M_PI/180.),
                        sin(angled*M_PI/180.),  cos(angled*M_PI/180.) };
  vector< double > gridu,
                   gridv;
  {
    // set bounding box min/max u/v
    double bbox[4] = { 1.e32, -1.e32, 1.e32, -1.e32 };

    for (unsigned i=0; i<(unsigned) points.size(); ++i) {
      const double u = (points[i].first-pc)^pu,
                   v = (points[i].first-pc)^pv,
                   w = (points[i].first-pc)^pn;
      points[i].second = PNT( r[0]*u+r[1]*v, r[2]*u+r[3]*v, w );
      bbox[0] = min(bbox[0],points[i].second[0]);
      bbox[1] = max(bbox[1],points[i].second[0]);
      bbox[2] = min(bbox[2],points[i].second[1]);
      bbox[3] = max(bbox[3],points[i].second[1]);
    }

    // set u/v grid
    gridu.assign(1,bbox[0]-du);
    gridv.assign(1,bbox[2]-dv);
    while (gridu.back()<bbox[1]+du)  gridu.push_back(gridu.back()+du);
    while (gridv.back()<bbox[3]+dv)  gridv.push_back(gridv.back()+dv);
  }
  cout << "info: project at rotation angle." << endl;


  cout << "info: split triangles..." << endl;
  // u coordinates splitting
  pair< unsigned,unsigned > splits(0,0);
  for (bool ok=false; !ok; ) {
    ok = true;
    const unsigned Ntri = (unsigned) tri.size();
    cout << "info: new pass, triangles to check: " << Ntri << "..." << endl;
    for (unsigned i=0; i<Ntri; ++i) {
      for (int e=0; e<3; ++e) {
        const int A = tri[i].v[NEXT(e)],  // edge node 1
                  B = tri[i].v[PREV(e)],  // edge node 2
                  C = tri[i].v[e],        // edge-opposite node
                  D = points.size();      // (possible) intersection node
        const double w = trisplit(gridu,points[A].second[0],points[B].second[0]);
        if (w>0. && w<1.) {
          points.push_back(pair< PNT,PNT >( points[A].first *(1.-w) + points[B].first *w,
                                            points[A].second*(1.-w) + points[B].second*w ));
          isedge.push_back(isedge[A] && isedge[B]);
          const double wacc = tri[i].w;
          tri.push_back(TRI());
          TRI& t1 = tri[i];
          TRI& t2 = tri.back();
          t1.v = Vec< 3,int >(C,A,D);  t1.i = tri[i].i;  t1.a = tri[i].a;  t1.w = wacc*w;
          t2.v = Vec< 3,int >(C,D,B);  t2.i = tri[i].i;  t2.a = tri[i].a;  t2.w = wacc*(1.-w);
          if (isedge.back()) {
            edge.erase(EDGE(A,B));
            edge.insert(EDGE(A,D));
            edge.insert(EDGE(B,D));
          }
          if (points[A].second[0]<points[B].second[0])  ++splits.first;
          else                                          ++splits.second;
          ok = false;
          break;
        }

      }  // for all edges
    }  // for all triangles
  }
  cout << "info: u splits: " << splits.first << '+' << splits.second << endl;

  // v coordinates splitting
  splits = pair< unsigned,unsigned >(0,0);
  for (bool ok=false; !ok; ) {
    ok = true;
    const unsigned Ntri = (unsigned) tri.size();
    cout << "info: new pass, triangles to check: " << Ntri << "..." << endl;
    for (unsigned i=0; i<Ntri; ++i) {
      for (int e=0; e<3; ++e) {
        const int A = tri[i].v[NEXT(e)],  // edge node 1
                  B = tri[i].v[PREV(e)],  // edge node 2
                  C = tri[i].v[e],        // edge-opposite node
                  D = points.size();      // (possible) intersection node
        const double w = trisplit(gridv,points[A].second[1],points[B].second[1]);
        if (w>0. && w<1.) {
          points.push_back(pair< PNT,PNT >( points[A].first *(1.-w) + points[B].first *w,
                                            points[A].second*(1.-w) + points[B].second*w ));
          isedge.push_back(isedge[A] && isedge[B]);
          const double wacc = tri[i].w;
          tri.push_back(TRI());
          TRI& t1 = tri[i];
          TRI& t2 = tri.back();
          t1.v = Vec< 3,int >(C,A,D);  t1.i = tri[i].i;  t1.a = tri[i].a;  t1.w = wacc*w;
          t2.v = Vec< 3,int >(C,D,B);  t2.i = tri[i].i;  t2.a = tri[i].a;  t2.w = wacc*(1.-w);
          if (isedge.back()) {
            edge.erase(EDGE(A,B));
            edge.insert(EDGE(A,D));
            edge.insert(EDGE(B,D));
          }
          if (points[A].second[1]<points[B].second[1])  ++splits.first;
          else                                          ++splits.second;
          ok = false;
          break;
        }

      }  // for all edges
    }  // for all triangles
  }
  cout << "info: v splits: " << splits.first << '+' << splits.second << endl;
  cout << "info: split triangles." << endl;


  cout << "info: build table..." << endl;
  vector< MNODE > table(gridu.size()*gridv.size());
  {
    // summary setup
    boost::progress_display pbar(table.size());
    unsigned nin = 0, nout = 0,
             nrcvy = 0, nfail = 0;
    vector< double > wsum(mold.e(0),0.);

    // build entry per quadrilateral
    for (unsigned i=0; i<(unsigned) gridu.size(); ++i) {
      for (unsigned j=0; j<(unsigned) gridv.size(); ++j, ++pbar) {
        MNODE& n = table[j*gridu.size() + i];

        // quadrilateral min/max x/y and center
        const double quad[4] = { gridu[i],gridu[i+1], gridv[j],gridv[j+1] };
        n.qcenter = PNT((quad[0]+quad[1])*.5,(quad[2]+quad[3])*.5,0.);
        n.nucsite_xyz = n.qcenter;

        // collect indices and weights withing quadrilateral...
        // ... and accumulate their contributions through a std::set< WPAIR >
        set< WPAIR,WPAIRLT > setwpair;
        vector< int > qedge;
        for (vector< TRI >::const_iterator t=tri.begin(); t!=tri.end(); ++t) {
          const PNT tc(( points[ t->v[0] ].second
                       + points[ t->v[1] ].second
                       + points[ t->v[2] ].second )/3.);
          if (tc[0]>quad[0] && tc[0]<quad[1] && tc[1]>quad[2] && tc[1]<quad[3]) {
            pair< set< WPAIR,WPAIRLT >::iterator,bool > r = setwpair.insert(WPAIR(t->i,t->w));
            if (!r.second /*if this index was already present*/) {
              const double wacc = (r.first)->second + t->w;
              setwpair.erase(r.first);
              setwpair.insert(WPAIR(t->i,wacc));
            }
            wsum[t->i] += t->w;
            if (isedge[(t->v)[0]])  qedge.push_back((t->v)[0]);
            if (isedge[(t->v)[1]])  qedge.push_back((t->v)[1]);
            if (isedge[(t->v)[2]])  qedge.push_back((t->v)[2]);
          }
        }
        nout += (setwpair.size()? 0:1);
        if (!setwpair.size())
          continue;
        for (set< WPAIR,WPAIRLT >::const_iterator p=setwpair.begin(); p!=setwpair.end(); ++p) {
          n.ni.push_back(p->first);
          n.nw.push_back(p->second);
        }

        // project onto xyz space by finding contaning triangle...
        // ... and interpolate using laplacian coordinates
        bool found = false;
        for (vector< TRI >::const_iterator t=tri.begin(); t!=tri.end() && !found; ++t) {
          const PNT &A2 = points[ t->v[0] ].second, &A3 = points[ t->v[0] ].first,
                    &B2 = points[ t->v[1] ].second, &B3 = points[ t->v[1] ].first,
                    &C2 = points[ t->v[2] ].second, &C3 = points[ t->v[2] ].first;
          if (tricheck(n.qcenter[0],A2[0],B2[0],C2[0],
                       n.qcenter[1],A2[1],B2[1],C2[1])) {
            const double a  = len(trinorm(A2,B2,C2)),
                         ka = len(trinorm(B2,C2,n.qcenter))/a,
                         kb = len(trinorm(C2,A2,n.qcenter))/a,
                         kc = len(trinorm(A2,B2,n.qcenter))/a;
            n.nucsite_xyz = A3*ka + B3*kb + C3*kc;
            n.nucsite_uv  = n.qcenter;
            found = true;
          }
        }
        nin  += (found? 1:0);
        nout += (found? 0:1);

        // if entry has contributions, but does not belong to any element...
        // ... shift to closest edge node inside quadrilateral
        if (!found) {
          sort(qedge.begin(),qedge.end());
          qedge.erase(unique(qedge.begin(),qedge.end()),qedge.end());
          double mind = sqrt(du*du+dv*dv)+1.;
          int    mini = -1;
          for (vector< int >::const_iterator i=qedge.begin(); i!=qedge.end(); ++i) {
            const double d = dist(n.nucsite_xyz,points[*i].second);
            if (d<mind) {
              mind = d;
              mini = *i;
            }
          }
          if (mini>=0) {
            n.nucsite_xyz = points[mini].first;
            const double u = (n.nucsite_xyz-pc)^pu,
                         v = (n.nucsite_xyz-pc)^pv,
                         w = (n.nucsite_xyz-pc)^pn;
            n.nucsite_uv = PNT( r[0]*u+r[1]*v, r[2]*u+r[3]*v, w );
          }
          else {
            n.ni.clear();
            n.nw.clear();
          }
          nin   += (mini<0? 0:1);  nout  += (mini<0? 1:0);
          nrcvy += (mini<0? 0:1);  nfail += (mini<0? 1:0);
        }

        // update (accumulate) quadrilateral area (intersected with zone)
        n.qarea = 0.;
        for (size_t j=0; j<n.ni.size(); ++j)
          n.qarea += tri[ n.ni[j] ].a * n.nw[j];

      }
    }

    // summary
    cout << "info: entries in/out (recovered/failed): " << nin << '/' << nout << " (" << nrcvy  << '/' << nfail << ')' << endl;
    cout << "info: element sum check, normalized [1.]: " << accumulate(wsum.begin(),wsum.end(),0.)/(double) wsum.size() << endl;

    // remove entries with no contribution
    for (int i=(int) table.size()-1; i>=0; --i)
      if (!table[i].ni.size())
        table.erase(table.begin()+i);
  }
  cout << "info: build table." << endl;


  {
    const string fn = "surfmap.table.xml";
    ofstream f;
    if (!ifstream(fn.c_str())) {
      f.open(fn.c_str(),ios::trunc);
      f << "<?xml version=\"1.0\" encoding=\"utf-8\"?>" << endl
        << "<!DOCTYPE ec [" << endl
        << " <!ENTITY year \"2009\">" << endl
        << " <!ENTITY author \"Vrije Universiteit Brussel\">" << endl
        << " <!ENTITY author-note \"&author;, &year;\">" << endl
        << "]>" << endl;
    }
    else {
      f.open(fn.c_str(),ios::app);
    }

    cout << "info: dump table (\"" << fn << "\")..." << endl;
    f << "<surfmapdb>" << endl
      << "  <surfmap zonename=\"" << zonename << "\" du=\"" << du << "\" dv=\"" << dv << "\" angled=\"" << angled << "\" n=\"" << table.size() << "\">" << endl;
    for (vector< MNODE >::const_iterator n=table.begin(); n!=table.end(); ++n) {
      f << "    <n"
        << " x=\"" << (n->nucsite_xyz)[0] << "\""
        << " y=\"" << (n->nucsite_xyz)[1] << "\""
        << " z=\"" << (n->nucsite_xyz)[2] << "\""
        << " e=\"" << (n->ni).size() << "\">" << endl;
      for (unsigned i=0; i<(n->ni).size(); ++i)
        f << "      <e i=\"" << (n->ni)[i] << "\" w=\"" << (n->nw)[i] << "\" />" << endl;
      f << "    </n>" << endl;
    }
    f << "  </surfmap>" << endl
      << "</surfmapdb>" << endl;
    cout << "info: dump table." << endl;
  }

  {
    const string fn2 = "surfmap.table.uv.plt";
    const string fn3 = "surfmap.table.xyz.plt";
    ofstream f2;
    ofstream f3;
    if (!ifstream(fn2.c_str()) || !ifstream(fn3.c_str())) {
      f2.open(fn2.c_str(),ios::trunc);
      f3.open(fn3.c_str(),ios::trunc);
      f2 << "VARIABLES = \"u\" \"v\"       \"weight\" \"area\"" << endl;
      f3 << "VARIABLES = \"x\" \"y\" \"z\" \"weight\" \"area\"" << endl;
    }
    cout << "info: dump table (\"" << fn2 << "\" and \"" << fn3 << "\")..." << endl;

    // triangulated (original) mesh
    for (unsigned i=0; i<(unsigned) table.size(); ++i) {
      const MNODE &N = table[i];

      f2 << "ZONE T=\"site " << i+1 << " contributions\""
         << " N=" << mold.n() << " E=" << mold.e(0)
         << " ZONETYPE=FETRIANGLE"
         << " DATAPACKING=BLOCK"
         << " VARLOCATION=([3-4]=CELLCENTERED)"
         << (i? " VARSHARELIST=([1-2,4]=1)":"") << endl;
      f3 << "ZONE T=\"site " << i+1 << " contributions\""
         << " N=" << mold.n() << " E=" << mold.e(0)
         << " ZONETYPE=FETRIANGLE"
         << " DATAPACKING=BLOCK"
         << " VARLOCATION=([4-5]=CELLCENTERED)"
         << (i? " VARSHARELIST=([1-3,5]=1)":"") << endl;

      for (unsigned d=0; d<2 && !i; ++d) {
        for (unsigned j=0; j<mold.n(); ++j)
          f2 << ' ' << points[j].second[d] << ((j+1)%500? ' ':'\n');
        f2 << endl;
      }
      for (unsigned d=0; d<3 && !i; ++d) {
        for (unsigned j=0; j<mold.n(); ++j)
          f3 << ' ' << points[j].first[d] << ((j+1)%500? ' ':'\n');
        f3 << endl;
      }
      for (int j=0; j<(int) mold.e(0); ++j) {
        double w = 0.;
        for (unsigned k=0; k<N.ni.size(); ++k)
          w += N.ni[k]==j? N.nw[k] : 0.;
        f2 << ' ' << w << ((j+1)%500? ' ':'\n');
        f3 << ' ' << w << ((j+1)%500? ' ':'\n');
      }
      for (unsigned j=0; j<mold.e(0) && !i; ++j) {
        f2 << ' ' << tri[j].a << ((j+1)%500? ' ':'\n');
        f3 << ' ' << tri[j].a << ((j+1)%500? ' ':'\n');
      }
      f2 << endl;
      f3 << endl;

      for (int j=0; j<(int) mold.e(0); ++j) {
        f2 << mold.vz[0].e2n[j].n[0]+1 << ' ' << mold.vz[0].e2n[j].n[1]+1 << ' ' << mold.vz[0].e2n[j].n[2]+1 << endl;
        f3 << mold.vz[0].e2n[j].n[0]+1 << ' ' << mold.vz[0].e2n[j].n[1]+1 << ' ' << mold.vz[0].e2n[j].n[2]+1 << endl;
      }
    }

    // nucleation site and size/quadrilaterals (contribution area)
    for (unsigned i=0; i<(unsigned) table.size(); ++i) {
      const MNODE &N = table[i];

      f2 << "ZONE T=\"site " << i+1 << "\" I=1"
         << " ZONETYPE=ORDERED"
         << " DATAPACKING=BLOCK"
         << " VARLOCATION=([4-5]=CELLCENTERED)" << endl
         << N.nucsite_uv[0] << endl
         << N.nucsite_uv[1] << endl
         << "0." << ' ' << N.qarea << endl;
      f3 << "ZONE T=\"site " << i+1 << "\" I=1"
         << " ZONETYPE=ORDERED"
         << " DATAPACKING=BLOCK"
         << " VARLOCATION=([3-4]=CELLCENTERED)" << endl
         << N.nucsite_xyz << endl
         << "0." << ' ' << N.qarea << endl;
      f2 << "ZONE T=\"quad " << i+1 << "\" N=4 E=1"
         << " ZONETYPE=FEQUADRILATERAL"
         << " DATAPACKING=BLOCK"
         << " VARLOCATION=([3-4]=CELLCENTERED)" << endl
         << N.qcenter[0]-du*.5 << ' ' << N.qcenter[0]-du*.5 << ' ' << N.qcenter[0]+du*.5 << ' ' << N.qcenter[0]+du*.5 << endl
         << N.qcenter[1]-dv*.5 << ' ' << N.qcenter[1]+dv*.5 << ' ' << N.qcenter[1]+dv*.5 << ' ' << N.qcenter[1]-dv*.5 << endl
         << "0." << ' ' << N.qarea << endl
         << "1 2 3 4" << endl;

    }

    cout << "info: dump table. (\"" << fn3 << "\")." << endl;
  }
}


double t_surfmap::trisplit(const vector< double >& vu, const double& u1, const double& u2)
{
  const double eps = 1.e-8;
  vector< double >::const_iterator lb = lower_bound(vu.begin(),vu.end(),min(u1,u2)+1.e-10);
  const double r = lb!=vu.end()? (*lb-min(u1,u2))/abs(u2-u1) : -1.;
  return (r>eps && r<1.-eps && u1<u2 && *lb<u2 && *lb>u1? r    :  // direct (u2-u1 is positive)
         (r>eps && r<1.-eps && u2<u1 && *lb<u1 && *lb>u2? 1.-r :  // reverse edge
                                                         -10. ));
}


bool t_surfmap::tricheck( const double& x, const double& xa, const double& xb, const double& xc,
                          const double& y, const double& ya, const double& yb, const double& yc )
{
  const Vec< 2,double > v0(xb-xa,yb-ya),
                        v1(xc-xa,yc-ya),
                        v2(x -xa,y -ya);
  const double dot00=v0^v0, dot01=v0^v1, dot02=v0^v2,
                            dot11=v1^v1, dot12=v1^v2;
  const double f = 1./(dot00*dot11 - dot01*dot01),
               u = (dot11*dot02 - dot01*dot12) * f,
               v = (dot00*dot12 - dot01*dot02) * f;
  return (u>=0) && (v>=0) && (u+v<=1);
}
