
#include <sstream>
#include "mfactory.h"
#include "t_simplex.h"

using namespace std;
using namespace m;


Register< mtransform,t_simplex > mt_simplex("-tzs","[str:...] split zone elements into simplex");


void t_simplex::transform(GetPot& o, mmesh& mold, const XMLNode& x)
{
  cout << "::simplex..." << endl;

  // find zones to simplex
  vector< unsigned > zs;
  string str = o.get(o.inc_cursor(),"");
  replace(str.begin(),str.end(),':',' ');
  istringstream is(str);
  string zname;
  while (is >> zname) {
    for (unsigned j=0; j<mold.vz.size(); ++j)
      if (mold.vz[j].n==zname) {
        zs.push_back(j);
        break;
      }
  }

  for (vector< unsigned >::iterator i=zs.begin(); i!=zs.end(); ++i) {
    mzone& z = mold.vz[*i];

    // generate simplexes for this zone (accumulator)
    vector< melem > vs;
    for (unsigned c=0; c<z.e2n.size(); ++c) {
      vector< melem > es = simplex(z.t,z.e2n[c].n);  // get element simplexes
      vs.insert(vs.end(),es.begin(),es.end());         // append to new connectivity
    }

    // setup new zone information
    if (vs.size()) {
      cout << "::simplex on zone \"" << z.n << "\" [e]: " << z.e2n.size() << '>' << vs.size() << endl;
      z.t = (z.d()==2? FETRIANGLE : (z.d()==3? FETETRAHEDRON : ORDERED ));
      z.e2n.swap(vs);
    }


  }
  cout << "::simplex." << endl;
}


vector< melem > t_simplex::simplex(const mtype& t, const vector< unsigned >& en)
{
  /*
   * this is an implementation of the paper:
   * "How to subdivide pyramids, prisms, and hexahedra into tetrahedra,"
   * J. Dompierre, P. Labbe, M.-G. Vallet, R. Camarero
   * Centre de recherche en calcul applique (CERCA)
   * 5160, boul. Decarie, bureau 400, Montreal, QC, H3X 2H9, Canada
   * [julien|paul|vallet|ricardo]@cerca.umontreal.ca
   * @misc{ dompierre99how,
   *        author = "J. Dompierre and P. Labbe and M. Vallet and R. Camarero",
   *        title = "How to subdivide pyramids, prisms and hexahedra into tetrahedra",
   *        text = "Dompierre J., Labbe P., Vallet M., Camarero R. How to subdivide pyramids,
   *          prisms and hexahedra into tetrahedra. Proceedings of the 8th International
   *          Meshing Roundtable, pp. 195--204. Sandia National Laboratories, Albuquerque,
   *          NM, 1999",
   *        year = "1999",
   *        url = "citeseer.ist.psu.edu/dompierre99how.html" }
   */
  melem simplex;
  simplex.n.assign(t==ORDERED?         1:  // no split
                  (t==FELINESEG?       3:  // no split
                  (t==FETRIANGLE?      3:  // no split
                  (t==FEQUADRILATERAL? 3:  // into triangle
                  (t==FETETRAHEDRON?   3:  // no split
                  (t==FEBRICK?         4:  // into tetrahedra
                  (t==FEPOLYGON?       3:  // no split (yet)
                  (t==FEPOLYHEDRON?    3:  // no split (yet)
                  (t==PRISM3?          4:  // into tetrahedra
                  (t==PYRAMID4?        4:  // into tetrahedra
                                       1 ))))))))),0);

  switch (t) {
    case FEQUADRILATERAL: {
      vector< melem > vsplit(2,simplex);

      // (gambit convection is the same)
      melem& s0 = vsplit[0];
      melem& s1 = vsplit[1];
      if (min(en[0],en[2])<min(en[1],en[3])) {
        s0.n[0] = en[0];  s0.n[1] = en[1];  s0.n[2] = en[2];
        s1.n[0] = en[0];  s1.n[1] = en[2];  s1.n[2] = en[3];
      }
      else {
        s0.n[0] = en[1];  s0.n[1] = en[2];  s0.n[2] = en[3];
        s1.n[0] = en[1];  s1.n[1] = en[3];  s1.n[2] = en[0];
      }

      return vsplit;
    } break;
    case FEBRICK: {
      vector< melem > vsplit(6,simplex);
      melem& s0 = vsplit[0];
      melem& s1 = vsplit[1];
      melem& s2 = vsplit[2];
      melem& s3 = vsplit[3];
      melem& s4 = vsplit[4];
      melem& s5 = vsplit[5];

      // nodal identifiers (by node numbering)
      vector< unsigned > V(8,0);

      // find minimum node identifier and number V's with indirect numbering (see paper)
      // (gambit convection is different, swapping nodes 2-3 and 6-7 (0-based))
#define n0 en[0]
#define n1 en[1]
#define n2 en[3]
#define n3 en[2]
#define n4 en[4]
#define n5 en[5]
#define n6 en[7]
#define n7 en[6]
      const unsigned M = *min_element(en.begin(),en.end());
      if      (n0==M) { V[0] = n0;  V[1] = n1;  V[2] = n2;  V[3] = n3;  V[4] = n4;  V[5] = n5;  V[6] = n6;  V[7] = n7; }
      else if (n1==M) { V[0] = n1;  V[1] = n0;  V[2] = n4;  V[3] = n5;  V[4] = n2;  V[5] = n3;  V[6] = n7;  V[7] = n6; }
      else if (n2==M) { V[0] = n2;  V[1] = n1;  V[2] = n5;  V[3] = n6;  V[4] = n3;  V[5] = n0;  V[6] = n4;  V[7] = n7; }
      else if (n3==M) { V[0] = n3;  V[1] = n0;  V[2] = n1;  V[3] = n2;  V[4] = n7;  V[5] = n4;  V[6] = n5;  V[7] = n6; }
      else if (n4==M) { V[0] = n4;  V[1] = n0;  V[2] = n3;  V[3] = n7;  V[4] = n5;  V[5] = n1;  V[6] = n2;  V[7] = n6; }
      else if (n5==M) { V[0] = n5;  V[1] = n1;  V[2] = n0;  V[3] = n4;  V[4] = n6;  V[5] = n2;  V[6] = n3;  V[7] = n7; }
      else if (n6==M) { V[0] = n6;  V[1] = n2;  V[2] = n1;  V[3] = n5;  V[4] = n7;  V[5] = n3;  V[6] = n0;  V[7] = n4; }
      else if (n7==M) { V[0] = n7;  V[1] = n3;  V[2] = n2;  V[3] = n6;  V[4] = n4;  V[5] = n0;  V[6] = n1;  V[7] = n5; }
#undef n0
#undef n1
#undef n2
#undef n3
#undef n4
#undef n5
#undef n6
#undef n7

      // set rotation bits, and rotate geometry 120/240 deg. by swapping
      unsigned s;
      const unsigned rb =
          (min(V[1],V[6])<min(V[2],V[5])? 4:0)   // for face 1
        + (min(V[3],V[6])<min(V[2],V[7])? 2:0)   // ... face 2
        + (min(V[4],V[6])<min(V[5],V[7])? 1:0);  // ... face 3
      if (rb==1 || rb==6) {
        s = V[1];  V[1] = V[4];  V[4] = V[3];  V[3] = s;
        s = V[5];  V[5] = V[7];  V[7] = V[2];  V[2] = s;
      }
      else if (rb==2 || rb==5) {
        s = V[1];  V[1] = V[3];  V[3] = V[4];  V[4] = s;
        s = V[5];  V[5] = V[2];  V[2] = V[7];  V[7] = s;
      }

      // split according to number of diagonal going through V[6]
      const unsigned N = (rb==0? 0 : (rb==7? 3 : (rb==1 || rb==2 || rb==4? 1:2)));
      if (N==0) {
        s0.n[0] = V[0];  s0.n[1] = V[1];  s0.n[2] = V[2];  s0.n[3] = V[5];
        s1.n[0] = V[0];  s1.n[1] = V[2];  s1.n[2] = V[7];  s1.n[3] = V[5];
        s2.n[0] = V[0];  s2.n[1] = V[2];  s2.n[2] = V[3];  s2.n[3] = V[7];
        s3.n[0] = V[0];  s3.n[1] = V[5];  s3.n[2] = V[7];  s3.n[3] = V[4];
        s4.n[0] = V[2];  s4.n[1] = V[7];  s4.n[2] = V[5];  s4.n[3] = V[6];
      }
      else if (N==1) {
        s0.n[0] = V[0];  s0.n[1] = V[5];  s0.n[2] = V[7];  s0.n[3] = V[4];
        s1.n[0] = V[0];  s1.n[1] = V[1];  s1.n[2] = V[7];  s1.n[3] = V[5];
        s2.n[0] = V[1];  s2.n[1] = V[6];  s2.n[2] = V[7];  s2.n[3] = V[5];
        s3.n[0] = V[0];  s3.n[1] = V[7];  s3.n[2] = V[2];  s3.n[3] = V[3];
        s4.n[0] = V[0];  s4.n[1] = V[7];  s4.n[2] = V[1];  s4.n[3] = V[2];
        s5.n[0] = V[1];  s5.n[1] = V[7];  s5.n[2] = V[6];  s5.n[3] = V[2];
      }
      else if (N==2) {
        s0.n[0] = V[0];  s0.n[1] = V[4];  s0.n[2] = V[5];  s0.n[3] = V[6];
        s1.n[0] = V[0];  s1.n[1] = V[3];  s1.n[2] = V[7];  s1.n[3] = V[6];
        s2.n[0] = V[0];  s2.n[1] = V[7];  s2.n[2] = V[4];  s2.n[3] = V[6];
        s3.n[0] = V[0];  s3.n[1] = V[1];  s3.n[2] = V[2];  s3.n[3] = V[5];
        s4.n[0] = V[0];  s4.n[1] = V[3];  s4.n[2] = V[6];  s4.n[3] = V[2];
        s5.n[0] = V[0];  s5.n[1] = V[6];  s5.n[2] = V[5];  s5.n[3] = V[2];
      }
      else if (N==3) {
        s0.n[0] = V[0];  s0.n[1] = V[2];  s0.n[2] = V[3];  s0.n[3] = V[6];
        s1.n[0] = V[0];  s1.n[1] = V[3];  s1.n[2] = V[7];  s1.n[3] = V[6];
        s2.n[0] = V[0];  s2.n[1] = V[7];  s2.n[2] = V[4];  s2.n[3] = V[6];
        s3.n[0] = V[0];  s3.n[1] = V[5];  s3.n[2] = V[6];  s3.n[3] = V[4];
        s4.n[0] = V[1];  s4.n[1] = V[5];  s4.n[2] = V[6];  s4.n[3] = V[0];
        s5.n[0] = V[1];  s5.n[1] = V[6];  s5.n[2] = V[2];  s5.n[3] = V[0];
      }

      if (!N)
        vsplit.pop_back();
      return vsplit;
    } break;
    case PRISM3: {
      vector< melem > vsplit(3,simplex);
      melem& s0 = vsplit[0];
      melem& s1 = vsplit[1];
      melem& s2 = vsplit[2];

      // nodal identifiers (by node numbering)
      vector< unsigned > V(6,0);

      // find minimum node identifier and number V's with indirect numbering (see paper)
      const unsigned M = *min_element(en.begin(),en.end());
      if      (en[0]==M) { V[0] = en[0];  V[1] = en[1];  V[2] = en[2];  V[3] = en[3];  V[4] = en[4];  V[5] = en[5]; }
      else if (en[1]==M) { V[0] = en[1];  V[1] = en[2];  V[2] = en[0];  V[3] = en[4];  V[4] = en[5];  V[5] = en[3]; }
      else if (en[2]==M) { V[0] = en[2];  V[1] = en[0];  V[2] = en[1];  V[3] = en[5];  V[4] = en[3];  V[5] = en[4]; }
      else if (en[3]==M) { V[0] = en[3];  V[1] = en[5];  V[2] = en[4];  V[3] = en[0];  V[4] = en[2];  V[5] = en[1]; }
      else if (en[4]==M) { V[0] = en[4];  V[1] = en[3];  V[2] = en[5];  V[3] = en[1];  V[4] = en[0];  V[5] = en[2]; }
      else if (en[5]==M) { V[0] = en[5];  V[1] = en[4];  V[2] = en[3];  V[3] = en[2];  V[4] = en[1];  V[5] = en[0]; }

      if (min(V[1],V[5])<min(V[2],V[4])) {
        s0.n[0] = V[0];  s0.n[1] = V[1];  s0.n[2] = V[2];  s0.n[3] = V[5];
        s1.n[0] = V[0];  s1.n[1] = V[1];  s1.n[2] = V[5];  s1.n[3] = V[4];
        s2.n[0] = V[0];  s2.n[1] = V[4];  s2.n[2] = V[5];  s2.n[3] = V[3];
      }
      else {
        s0.n[0] = V[0];  s0.n[1] = V[1];  s0.n[2] = V[2];  s0.n[3] = V[4];
        s1.n[0] = V[0];  s1.n[1] = V[4];  s1.n[2] = V[2];  s1.n[3] = V[5];
        s2.n[0] = V[0];  s2.n[1] = V[4];  s2.n[2] = V[5];  s2.n[3] = V[3];
      }

      return vsplit;
    } break;
    case PYRAMID4: {
      vector< melem > vsplit(2,simplex);
      melem& s0 = vsplit[0];
      melem& s1 = vsplit[1];

      // (gambit convection is different, swapping nodes 2-3 (0-based))
#define n0 en[0]
#define n1 en[1]
#define n2 en[3]
#define n3 en[2]
#define n4 en[4]
      if (min(n0,n2)<min(n1,n3)) {
        s0.n[0] = n0;  s0.n[1] = n1;  s0.n[2] = n2;  s0.n[3] = n4;
        s1.n[0] = n0;  s1.n[1] = n2;  s1.n[2] = n3;  s1.n[3] = n4;
      }
      else {
        s0.n[0] = n1;  s0.n[1] = n2;  s0.n[2] = n3;  s0.n[3] = n4;
        s1.n[0] = n1;  s1.n[1] = n3;  s1.n[2] = n0;  s1.n[3] = n4;
      }
#undef n0
#undef n1
#undef n2
#undef n3
#undef n4

      return vsplit;
    } break;
    case ORDERED:
    case FELINESEG:
    case FETRIANGLE:
    case FETETRAHEDRON:
    case FEPOLYGON:
    case FEPOLYHEDRON:
    default: break;
  }

  return vector< melem >();
}

