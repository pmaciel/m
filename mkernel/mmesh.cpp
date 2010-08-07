
#include <iostream>
#include <algorithm>
#include "mmesh.h"

using namespace std;


namespace m {


unsigned mmesh::d() const
{
  string trail(3,'?');  // last characters of first 3 names
  for (unsigned i=0; i<3 && i<vn.size(); ++i)
    if (vn[i].length())
      trail[i] = *(vn[i].end()-1);

  unsigned d = 0;  // dimension
  if (trail[0]!='?')                 d=1;
  if (d==1 && trail[1]==trail[0]+1)  d=2;
  if (d==2 && trail[2]==trail[0]+2)  d=3;
  return d;

  /*
  works well, but might be difficult to understand
  if (trail[0]!='?') { if (trail[1]==trail[0]+1) { if (trail[2]==trail[0]+2) {
  return 3; } else return 2; } else return 1; } else return 0;
  */
}


vector< bool > mmesh::vvectors() const
{
  // last characters of variables names, with 2 extra to simplify check
  string trail(vn.size()+2,'?');
  for (unsigned i=0; i<vn.size(); ++i)
    if (vn[i].length())
      trail[i] = *(vn[i].end()-1);

  // check names
  vector< bool > r(vn.size(),false);
  const unsigned dim = d();
  if (dim<2)
    return r;

  for (unsigned i=0; i<r.size(); ++i) {
    r[i] = (trail[i]!='?');
    if (dim>1 && r[i])  r[i] = (trail[i+1]==trail[i]+1);
    if (dim>2 && r[i])  r[i] = (trail[i+2]==trail[i]+2);
    if (r[i]) {
      cout << "info: detected vector field, variables ";
      for (unsigned j=0; j<dim; ++j)
        cout << char('x'+j) << ":" << vn[i+j] << " ";
      cout << endl;
    }
  }
  return r;
}


void mmesh::merge(const mmesh& another)
{
/*
  // -- merge dimensions
  // -- merge coordinates, other variables are discarded --
  // -- merge nodes (renumber, removing duplicates) --
  // -- merge nodes (renumber element > nodes lists) --
  // -- create inverse connectivity --
  // -- merge elements (see only elements with renumbered nodes) --
*/
}


void mmesh::compress()
{
  // get used zones node indices and remove duplicates
  // B2A: new-to-old numbering
  // A2B: old-to-new numbering
  vector< unsigned > B2A,
                     A2B;
  for (vector< mzone >::const_iterator z=vz.begin(); z!=vz.end(); ++z)
    for (unsigned i=0; i<z->e2n.size(); ++i)
      B2A.insert(B2A.end(),
        z->e2n[i].n.begin(),z->e2n[i].n.end() );
  sort(B2A.begin(),B2A.end());
  B2A.erase(unique(B2A.begin(),B2A.end()),B2A.end());
  A2B.assign(n(),0);
  for (unsigned i=0; i<B2A.size(); ++i)
    A2B[ B2A[ i ] ] = i;

  // re-index point cloud values...
  for (unsigned i=0; i<vv.size(); ++i) {
    vector< double > vvv(B2A.size(),0.);
    for (unsigned j=0; j<B2A.size(); ++j)
      vvv[j] = vv[i][B2A[j]];
    vv[i].swap(vvv);
  }

  // and renumber connectivities
  for (vector< mzone >::iterator z=vz.begin(); z!=vz.end(); ++z) {
    for (unsigned i=0; i<z->e2n.size(); ++i)
      for (vector< unsigned >::iterator j=z->e2n[i].n.begin(); j!=z->e2n[i].n.end(); ++j)
        *j = A2B[*j];
  }
}


mmesh mmesh::extract(const string& zn)
{
  mmesh r;
  for (vector< mzone >::iterator z=vz.begin(); z!=vz.end(); ++z) {
    if (z->n==zn) {
      // get used zone node indices and remove duplicates
      // B2A: new-to-old numbering
      // A2B: old-to-new numbering
      vector< unsigned > B2A,
                         A2B;
      for (unsigned i=0; i<z->e2n.size(); ++i)
        B2A.insert(B2A.end(),
          z->e2n[i].n.begin(),z->e2n[i].n.end() );
      sort(B2A.begin(),B2A.end());
      B2A.erase(unique(B2A.begin(),B2A.end()),B2A.end());
      A2B.assign(vv[0].size(),0);
      for (unsigned i=0; i<B2A.size(); ++i)
        A2B[ B2A[ i ] ] = i;

      // re-index point cloud values...
      const unsigned Nnode = B2A.size();
      r.vn = vn;
      r.vv.assign(v(),vector< double >(Nnode));
      for (unsigned i=0; i<Nnode; ++i)
        for (unsigned j=0; j<v(); ++j)
          r.vv[j][i] = vv[j][B2A[i]];

      // and renumber (this zone) connectivity
      r.vz.assign(1,*z);   // copy given zone
      for (unsigned i=0; i<r.e(0); ++i)
        for (vector< unsigned >::iterator j=r.vz[0].e2n[i].n.begin(); j!=r.vz[0].e2n[i].n.end(); ++j)
          *j = A2B[*j];

      return r;
    }
  }

  cerr << "error: asking zone \"" << zn << "\", does not exist" << endl;
  throw 42;
  return r;
}


}  // namespace m

