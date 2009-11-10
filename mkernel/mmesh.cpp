
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
  vector< unsigned > B2A;
  vector< int >      A2B;
  for (vector< mzone >::const_iterator z=vz.begin(); z!=vz.end(); ++z)
    for (unsigned i=0; i<z->e2n.size(); ++i)
      B2A.insert(B2A.end(),
        z->e2n[i].n.begin(),z->e2n[i].n.end() );
  sort(B2A.begin(),B2A.end());
  B2A.erase(unique(B2A.begin(),B2A.end()),B2A.end());
  A2B.assign(n(),0);
  for (unsigned i=0; i<B2A.size(); ++i)
    A2B[ B2A[ i ] ] = (int) i;

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
        *j = (unsigned) A2B[*j];
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

      // point cloud names and values
      const unsigned Nnode = B2A.size();
      r.vn = vn;
      r.vv.assign(v(),vector< double >(Nnode));
      for (unsigned i=0; i<Nnode; ++i)
        for (unsigned j=0; j<v(); ++j)
          r.vv[j][i] = vv[j][B2A[i]];

      // set just this zone
      r.vz.assign(1,*z);   // copy given zone
      for (unsigned i=0; i<r.e(0); ++i)
        for (vector< unsigned >::iterator j=r.vz[0].e2n[i].n.begin(); j!=r.vz[0].e2n[i].n.end(); ++j)
          *j = A2B[*j];   // and renumber its connectivity

      return r;
    }
  }

  cerr << "error: asking zone \"" << zn << "\", does not exist" << endl;
  throw 42;
  return r;
}


mmesh mmesh::edge(const string& zn)
{
  if (z()<2) {
    cerr << "error: asking edge of \"" << zn << "\", but (z==" << z() << ")>=2" << endl;
    throw 42;
  }

  mmesh r;
  for (unsigned i=0; i<z(); ++i) {
    if (vz[i].n==zn) {

      /*
       * detect edge nodes by marking if they are used by our zone, and by
       * another one (with same dimensionality) as well
       */
      vector< bool > ninedge(n(),false);
      {
        vector< vector< unsigned > > ninzone(n());
        for (unsigned j=0; j<z(); ++j) {
          if (d(j)==d(i))
            for (vector< melem >::const_iterator e=vz[j].e2n.begin(); e!=vz[j].e2n.end(); ++e) {
              for (vector< unsigned >::const_iterator n=(*e).n.begin(); n!=(*e).n.end(); ++n)
                ninzone[ *n ].push_back( j );
            }
        }
        for (unsigned j=0; j<n(); ++j) {
          sort(ninzone[j].begin(),ninzone[j].end());
          ninzone[j].erase(unique(ninzone[j].begin(),ninzone[j].end()),ninzone[j].end());
          ninedge[j] = (ninzone[j].size()>=2 && count(ninzone[j].begin(),ninzone[j].end(),i));
        }
      }

      /*
       * list all elements using two or more of those nodes, that is, with one
       * side on the edge; the elements are sorted by zone
       */
      vector< vector< melem > > eelems(z());
      for (unsigned j=0; j<z(); ++j) {
        if (d(j)==d(i))
          for (vector< melem >::const_iterator e=vz[j].e2n.begin(); e!=vz[j].e2n.end(); ++e) {
            unsigned c = 0;
            for (vector< unsigned >::const_iterator n=(*e).n.begin(); n!=(*e).n.end(); ++n)
              c += ( ninedge[ *n ]? 1:0);
            if (c>1)
              eelems[j].push_back( *e );
          }
      }

      /*
       * create line segments by, for all edge elements of our zone, match it
       * against all other edge elements and keep their common sides (thus, they
       * share two nodes) - this way, only zero or one shared sides are possible
       */
      r.vz.resize(1);
      r.vz[0].n = vz[i].n + "_edge";
      for (vector< melem >::const_iterator ei=eelems[i].begin(); ei!=eelems[i].end(); ++ei) {
        for (unsigned j=0; j<z(); ++j)
          if (i!=j)
            for (vector< melem >::const_iterator ej=eelems[j].begin(); ej!=eelems[j].end(); ++ej) {
              melem seg;
              for (vector< unsigned >::const_iterator k=ei->n.begin(); k!=ei->n.end(); ++k)
                if (count(ej->n.begin(),ej->n.end(),*k))
                  seg.n.push_back(*k);
              if (seg.n.size()==2)
                r.vz[0].e2n.push_back(seg);
            }
      }
      if (!r.e(0)) {
        cerr << "error: asking edge of \"" << zn << "\", but no edge was found" << endl;
        throw 42;
      }
      cout << "info: detected edge elements: " << r.e(0) << endl;


      /*
       * sort line segments by begining with the first and loop around, locating
       * the next segment with the same node; also, line segments "go forward"
       */
      {
        const vector< melem > old = r.vz[0].e2n;
        for (vector< melem >::iterator e1=r.vz[0].e2n.begin(); e1!=r.vz[0].e2n.end()-1; ) {
          const unsigned a1 = (*e1).n[0];
          const unsigned b1 = (*e1).n[1];
          bool found = false;
          for (vector< melem >::const_iterator e2=old.begin(); e2!=old.end() && !found; ++e2) {
            const unsigned a2 = e2->n[0];
            const unsigned b2 = e2->n[1];
            if      (b1==a2 && a1!=b2) { found=true; (++e1)->n[0]=a2; e1->n[1]=b2; }
            else if (b1==b2 && a1!=a2) { found=true; (++e1)->n[0]=b2; e1->n[1]=a2; }
          }
          if (!found) {
            cerr << "error: edge doesn't form a closed loop" << endl;
            throw 42;
          }
        }
      }
      cout << "info: sorted edge elements: " << r.e(0) << endl;

      /*
       * copy point cloud, and compress mesh
       */
      r.vv = vv;
      r.vn = vn;
      r.compress();

      return r;
    }
  }

  cerr << "error: asking zone \"" << zn << "\", does not exist" << endl;
  throw 42;
  return r;
}


}  // namespace m

