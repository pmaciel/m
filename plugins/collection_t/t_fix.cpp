
#include <numeric>
#include "mfactory.h"
#include "t_fix.h"  // this class implementation
#include "../io_smurf/f_plt.h"  // use this class to write tecplot files
#include "cool/Element.hh"  // elements properties

using namespace std;
using namespace m;


Register< mtransform,t_fix > mt_fix("-tfix","very heavy super mesh fix");


void t_fix::transform(GetPot& o, mmesh& m, const XMLNode& x)
{
  const unsigned d = m.d();
  for (vector< mzone >::const_iterator iz=m.vz.begin(); iz!=m.vz.end(); ++iz) {
    if (d==2) {
      if      (iz->t==FELINESEG)      check< 2,FELINESEG     >(*iz,m.vn,m.vv);
      else if (iz->t==FETRIANGLE)     check< 2,FETRIANGLE    >(*iz,m.vn,m.vv);
    }
    else if (d==3) {
      if      (iz->t==FELINESEG)      check< 3,FELINESEG     >(*iz,m.vn,m.vv);
      else if (iz->t==FETRIANGLE)     check< 3,FETRIANGLE    >(*iz,m.vn,m.vv);
      else if (iz->t==FETETRAHEDRON)  check< 3,FETETRAHEDRON >(*iz,m.vn,m.vv);
    }
  }
}


/*
 * the templated method implementation is here because no code includes
 * t_fix.h (not even main program, because we have self-registration).
 * it is incorrect to use "using" in a header... and i want it :)
 */
template< int D, int T >
void t_fix::check(const mzone& z, const vector< string >& vn, const vector< vector< double > >& vv)
{
  using namespace ::COOLFluiD::Muffin;
  using ::COOLFluiD::Framework::Node;
  using ::COOLFluiD::Framework::Node;
  using ::COOLFluiD::Framework::State;

  // number of nodes, elements and states
  const unsigned I = vv.size()? vv[0].size():0;
  const unsigned J = z.e2n.size();
  const unsigned V = (unsigned) max< int >(0,vv.size()-D);

  // states derivatives (possibly), their names and nodal size (dual)
  vector< vector< double > > dstatedx,  dstatedy,  dstatedz;
  vector< CFreal >           ddx,       ddy,       ddz;
  if (V && D>0) { dstatedx.assign(V,vector< double >(I,0.));  ddx.assign(V,0.); }
  if (V && D>1) { dstatedy.assign(V,vector< double >(I,0.));  ddy.assign(V,0.); }
  if (V && D>2) { dstatedz.assign(V,vector< double >(I,0.));  ddx.assign(V,0.); }
  vector< string > newvn(vn);
  newvn.push_back("nodal_size");
  vector< double > nodal_size(I,0.);
  if (V && D>0) for (unsigned s=0; s<V; ++s) newvn.push_back(vn[D+s]+"_ddx");
  if (V && D>1) for (unsigned s=0; s<V; ++s) newvn.push_back(vn[D+s]+"_ddy");
  if (V && D>2) for (unsigned s=0; s<V; ++s) newvn.push_back(vn[D+s]+"_ddz");

  // allocate memory for element, Node*s and State*s
  AElement* e(z.t==FELINESEG?       (AElement*) new ElementLineseg(D)     :
             (z.t==FETRIANGLE?      (AElement*) new ElementTriangle(D)    :
             (z.t==FETETRAHEDRON?   (AElement*) new ElementTetrahedron(D) :
                                    (AElement*) NULL )));
  if (e==(AElement*) NULL) {
    cout << "Warning: zone \"" << z.n << "\" will not be checked because its type is not supported" << endl;
    return;
  }
  vector< Node* > enodes;
  for (CFuint i=0; i<e->N; ++i)
    enodes.push_back(new Node(0.,D));
  if (V) {
    for (CFuint i=0; i<e->S; ++i) {
      e->states[i].resize(V);
      e->states[i] = 0.;
    }
  }

  double ssum = 0.;
  double smin = 1.e99;
  double smax = 0.;
  for (vector< melem >::const_iterator ie=z.e2n.begin(); ie!=z.e2n.end(); ++ie) {

    // get element's nodes and states
    for (unsigned i=0; i<e->N; ++i)
      for (unsigned j=0; j<D; ++j)
        (*enodes[i])[j] = vv[j][ie->n[i]];


    // get element's states
    if (V) {
      for (unsigned i=0; i<e->N; ++i)
        for (unsigned j=0; j<V; ++j)
          e->states[i][j] = vv[D+j][ie->n[i]];
    }


    // build element
    e->element(enodes);


    // update size statistics
    ssum += e->s;
    smin = min(smin,e->s);
    smax = max(smax,e->s);


    // calculate and distribute derivatives weighted with nodal size
    const double w = (e->s/e->N);
    for (unsigned i=0; i<ie->n.size(); ++i)
      nodal_size[ie->n[i]] += w;
    for (unsigned i=0; i<e->N; ++i) {
      for (unsigned s=0; s<V; ++s) if (D>0) { ddx = e->dd(XX);  dstatedx[s][ie->n[i]] += ddx[s]*w; }
      for (unsigned s=0; s<V; ++s) if (D>1) { ddy = e->dd(YY);  dstatedy[s][ie->n[i]] += ddy[s]*w; }
      for (unsigned s=0; s<V; ++s) if (D>2) { ddz = e->dd(ZZ);  dstatedz[s][ie->n[i]] += ddz[s]*w; }
    }

  }


  // normalize derivatives with nodal size
  for (unsigned n=0; V && n<I; ++n) if (nodal_size[n]>1.e-20) for (unsigned s=0; s<V; ++s) {
        if (D>0)  dstatedx[s][n] /= nodal_size[n];
        if (D>1)  dstatedy[s][n] /= nodal_size[n];
        if (D>2)  dstatedz[s][n] /= nodal_size[n];
  }
  cout << "zone \"" << z.n << "\" d=" << z.d() << ": " << endl
       << "  size (sum e.): "     << ssum << endl
       << "  size (sum n.): "     << accumulate(nodal_size.begin(),nodal_size.end(),0.) << endl
       << "  element min. size: " << smin << endl
       << "  element max. size: " << smax << endl
       << "  element avg. size: " << (J? ssum/J:0.) << endl;


  // write a file per zone
  const string fn(string("fix_" + z.n + ".plt").c_str());
  ofstream f(fn.c_str());
  if (!f) {
    cerr << "error accessing file: \"" << fn << "\"" << endl;
    throw 42;
  }
  f.precision(15);

  TecZone tz(z.n, e->N, I, J, z.t, true /*isblock*/, false /*isshared*/);
  f << f_plt::setVariables(newvn) << endl
    << f_plt::setZoneHeader(tz) << endl;
  f_plt::writeZoneNodeValues(f,vv,tz.isblock);
  f_plt::writeZoneNodeValues(f,nodal_size);
  if (V && D>0)  f_plt::writeZoneNodeValues(f,dstatedx,tz.isblock);
  if (V && D>1)  f_plt::writeZoneNodeValues(f,dstatedy,tz.isblock);
  if (V && D>2)  f_plt::writeZoneNodeValues(f,dstatedz,tz.isblock);
  f_plt::writeZoneConnectivity(f,z.e2n,z.t);
  f.close();


  // deallocate memory ef element and Node*s
  delete e;
  for (CFuint i=0; i<enodes.size(); ++i)
    delete enodes[i];
}


