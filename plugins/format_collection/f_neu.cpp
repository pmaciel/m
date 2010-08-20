
#include "mfactory.h"
#include "f_neu.h"
#include <sstream>

using namespace std;
using namespace m;


Register< mfinput,f_neu > mf_neu(".neu","Fluent Gambit neutral input format");


// auxiliary type for gambit elements
struct GElement {
  vector< unsigned > n;  // nodes
  mtype t;               // type
  unsigned g;            // group
};


// auxiliary type for gambit groups
struct GGroup {
  string n;              // name
  vector< unsigned > e;  // elements
};


void f_neu::read(GetPot& o,mmesh& m)
{
  const string fn(o.get(o.inc_cursor(),""));
  ifstream f(fn.c_str());
  if (!f) {
    cerr << "error accessing file: \"" << fn << "\"" << endl;
    throw 42;
  }


  int dummy;
  unsigned NUMNP,
           NELEM,
           NGRPS,
           NBSETS,
           NDFCD /* not used*/,
           NDFVL;
  string s;

  // skip 6 initial lines
  for (int i=0; i<6; ++i)
    getline(f,s);

  // read number of points, elements, groups, sets and dimensions
  getline(f,s);
  istringstream ss(s);
  ss >> NUMNP >> NELEM >> NGRPS >> NBSETS >> NDFCD >> NDFVL;

  // set variable names (just coordinates)
  m.vn.resize(NDFVL);
  for (unsigned d=0; d<NDFVL; ++d)
    m.vn[d] = string(1,char('x'+d));

  // skip next two lines
  for (unsigned i=0; i<2; ++i)
    getline(f,s);

  // allocate and read coordinates
  m.vv.assign(NDFVL,vector< double >(NUMNP,0.));
  for (unsigned i=0; i<NUMNP; ++i) {
    getline(f,s);
    istringstream ss(s);
    unsigned ND;
    ss >> ND;
    ND--;
    for (unsigned d=0; d<NDFVL; ++d)
      ss >> m.vv[d][ND];
  }

  // skip next two lines
  for (unsigned i=0; i<2; ++i)
    getline(f,s);

  // allocate and read elements connectivity (structure is 1-based)
  vector< GElement > vgelems(NELEM);
  for (vector< GElement >::iterator e=vgelems.begin(); e!=vgelems.end(); ++e) {

    // element description
    unsigned NE, NTYPE, NDP, n;
    f >> NE >> NTYPE >> NDP;
    if      (NTYPE==2 && NDP==4) e->t=FEQUADRILATERAL;  // quadrilateral
    else if (NTYPE==3 && NDP==3) e->t=FETRIANGLE;       // triangle
    else if (NTYPE==4 && NDP==8) e->t=FEBRICK;          // brick
    else if (NTYPE==5 && NDP==6) e->t=PRISM3;           // wedge (prism)
    else if (NTYPE==6 && NDP==4) e->t=FETETRAHEDRON;    // tetrahedron
    else if (NTYPE==7 && NDP==5) e->t=PYRAMID4;         // pyramid
    else {
      cerr << "error: no support for element type/nodes " << NTYPE << "/" << NDP << endl;
      throw 42;
    }

    // get element nodes
    e->n.assign(NDP,0);
    for (unsigned j=0; j<NDP; ++j) {
      f >> n;
      e->n[j] = n-1;  // 0-based index
    }

    // finish the line
    getline(f,s);
  }

  getline(f,s);  // ENDOFSECTION

  // read element groups (volume)
  vector< GGroup > vggroups(NGRPS);
  vector< vector< unsigned > > nbelems(NGRPS,vector< unsigned >(PYRAMID4+1,0)); // nb. elements in group of type
  for (unsigned g=0; g<NGRPS; ++g) {
    string ELMMAT;
    unsigned NGP, NELGP, MTYP, NFLAGS, I;
    getline(f,s);  // ELEMENT GROUP...

    f >> s >> NGP >> s >> NELGP >> s >> MTYP >> s >> NFLAGS >> ELMMAT;
    vggroups[g].n = ELMMAT;
    vggroups[g].e.resize(NELGP);

    for (unsigned i=0; i<NFLAGS; ++i)
      f >> dummy;
    for (unsigned i=0; i<NELGP; ++i) {
      f >> I;
      vggroups[g].e[i] = I-1;           // set element index
      ++nbelems[g][ vgelems[ I-1 ].t ]; // update counter
    }

    getline(f,s);  // finish the line (read new line)
    getline(f,s);  // ENDOFSECTION
  }

  // assign elements to volume zones, distinguishing between element types
  m.vz.clear();
  for (unsigned g=0; g<NGRPS; ++g) {

    // count element types in group
    vector< mtype    > gelemstypes;
    vector< unsigned > gelemsnb;
    for (unsigned t=0; t<nbelems[g].size(); ++t)
      if (nbelems[g][t]>0) {
        gelemstypes.push_back((mtype) t);
        gelemsnb.push_back(nbelems[g][t]);
      }

    // add new zones
    for (unsigned i=0; i<gelemstypes.size(); ++i) {
      m.vz.push_back(mzone());
      mzone& z = m.vz.back();

      // set new zone name
      ostringstream name;
      name << vggroups[g].n;
      if (gelemstypes.size()>1)
        name << "_g" << g+1 << "_t" << gelemstypes[i];
      z.n = name.str();

      // set new zone type and dimensionality
      z.t = gelemstypes[i];
      if (!z.d()) {
        cerr << "error: error detecting zone dimensionality" << endl;
        throw 42;
      }

      // set new zone element-node connectivity
      z.e2n.resize(gelemsnb[i]);
      for (unsigned j=0, k=0; j<vggroups[g].e.size(); ++j)
        if (vgelems[ vggroups[g].e[j] ].t==gelemstypes[i])
          z.e2n[k++].n = vgelems[ vggroups[g].e[j] ].n;
    }
  }

  // read each boundary separately
  for (unsigned t=0; t<NBSETS; ++t) {

    string NAME;
    int ITYPE, NENTRY, NVALUES, IBCODE1, IBCODE2, IBCODE3, IBCODE4, IBCODE5;

    // read header
    getline(f,s);  // BOUNDARY CONDITIONS...
    getline(f,s);  // header
    istringstream ss(s);
    ss >> NAME >> ITYPE >> NENTRY >> NVALUES >> IBCODE1 >> IBCODE2 >> IBCODE3 >> IBCODE4 >> IBCODE5;
    if (ITYPE!=1) {
      cerr << "error: supports only boundary condition data 1 (element/cell): page C-11 of user's guide" << endl;
      throw 42;
    }
    if (IBCODE1!=6) {
      cerr << "error: supports only IBCODE1 6 (ELEMENT_SIDE)" << endl;
      throw 42;
    }

    // boundary connectivity here
    vector< GElement > e2n(NENTRY);

    // read boundary elements connectivity
    vector< unsigned > nbelems(PYRAMID4+1,0); // nb. elements per type
    for (int i=0; i<NENTRY; ++i) {
      int ELEM, ETYPE, FACE;
      f >> ELEM >> ETYPE >> FACE;

      // element nodes and face nodes/type
      const vector< unsigned >& en = vgelems[ ELEM-1 ].n;
      vector< unsigned >& fn = e2n[i].n;
      mtype&              ft = e2n[i].t;
      if      (ETYPE==2 && en.size()==4)  ft = FELINESEG;        // quadrilateral faces
      else if (ETYPE==3 && en.size()==3)  ft = FELINESEG;        // triangle faces
      else if (ETYPE==4 && en.size()==8)  ft = FEQUADRILATERAL;  // brick faces
      else if (ETYPE==5 && en.size()==6)  ft = (FACE<4? FEQUADRILATERAL : FETRIANGLE);  // wedge (prism) faces
      else if (ETYPE==6 && en.size()==4)  ft = FETRIANGLE;       // tetrahedron faces
      else if (ETYPE==7 && en.size()==5)  ft = (FACE<2? FEQUADRILATERAL : FETRIANGLE);  // pyramid faces
      else {
        cerr << "error: reference for an unexpected volume element" << endl;
        throw 42;
      }
      ++nbelems[(unsigned) ft];  // add a face element of this type
      fn.assign(ft==FELINESEG?       2:
               (ft==FEQUADRILATERAL? 4:
               (ft==FETRIANGLE?      3:
                                     1 )), 0);

      if (ETYPE==2) { // quadrilateral faces
        switch (FACE) {
          case 1: fn[0] = en[0]; fn[1] = en[1]; break;
          case 2: fn[0] = en[1]; fn[1] = en[2]; break;
          case 3: fn[0] = en[2]; fn[1] = en[3]; break;
          case 4: fn[0] = en[3]; fn[1] = en[0]; break;
        }
      }
      else if (ETYPE==3) { // triangle faces
        switch (FACE) {
          case 1: fn[0] = en[0]; fn[1] = en[1]; break;
          case 2: fn[0] = en[1]; fn[1] = en[2]; break;
          case 3: fn[0] = en[2]; fn[1] = en[0]; break;
        }
      }
      else if (ETYPE==4) { // brick faces
        switch (FACE) {
          case 1: fn[0] = en[0]; fn[1] = en[1]; fn[2] = en[5]; fn[3] = en[4]; break;
          case 2: fn[0] = en[1]; fn[1] = en[3]; fn[2] = en[7]; fn[3] = en[5]; break;
          case 3: fn[0] = en[3]; fn[1] = en[2]; fn[2] = en[6]; fn[3] = en[7]; break;
          case 4: fn[0] = en[2]; fn[1] = en[0]; fn[2] = en[4]; fn[3] = en[6]; break;
          case 5: fn[0] = en[1]; fn[1] = en[0]; fn[2] = en[2]; fn[3] = en[3]; break;
          case 6: fn[0] = en[4]; fn[1] = en[5]; fn[2] = en[7]; fn[3] = en[6]; break;
        }
      }
      else if (ETYPE==5) { // wedge (prism) faces
        switch (FACE) {
          case 1: fn[0] = en[0]; fn[1] = en[1]; fn[2] = en[4]; fn[3] = en[3]; break;
          case 2: fn[0] = en[1]; fn[1] = en[2]; fn[2] = en[5]; fn[3] = en[4]; break;
          case 3: fn[0] = en[2]; fn[1] = en[0]; fn[2] = en[3]; fn[3] = en[5]; break;
          case 4: fn[0] = en[0]; fn[1] = en[2]; fn[2] = en[1]; break;
          case 5: fn[0] = en[3]; fn[1] = en[4]; fn[2] = en[5]; break;
        }
      }
      else if (ETYPE==6) { // tetrahedron faces
        switch (FACE) {
          case 1: fn[0] = en[1]; fn[1] = en[0]; fn[2] = en[2]; break;
          case 2: fn[0] = en[0]; fn[1] = en[1]; fn[2] = en[3]; break;
          case 3: fn[0] = en[1]; fn[1] = en[2]; fn[2] = en[3]; break;
          case 4: fn[0] = en[2]; fn[1] = en[0]; fn[2] = en[3]; break;
        }
      }
      else if (ETYPE==7) { // pyramid faces
        switch (FACE) {
          case 1: fn[0] = en[0]; fn[1] = en[2]; fn[2] = en[3]; fn[3] = en[1]; break;
          case 2: fn[0] = en[0]; fn[1] = en[1]; fn[2] = en[4]; break;
          case 3: fn[0] = en[1]; fn[1] = en[3]; fn[2] = en[4]; break;
          case 4: fn[0] = en[3]; fn[1] = en[2]; fn[2] = en[4]; break;
          case 5: fn[0] = en[2]; fn[1] = en[0]; fn[2] = en[4]; break;
        }
      }
      getline(f,s);  // finish the line (read new line)

    }
    getline(f,s);  // ENDOFSECTION


    // add a new zone, splitting according to type (if necessary)
    // count element types
    vector< mtype    > felemstypes;
    vector< unsigned > felemsnb;
    for (unsigned t=0; t<nbelems.size(); ++t)
      if (nbelems[t]>0) {
        felemstypes.push_back((mtype) t);
        felemsnb.push_back(nbelems[t]);
      }

    // set new zones, distinguishing different element types
    for (unsigned i=0; i<felemstypes.size(); ++i) {
      m.vz.push_back(mzone());
      mzone& z = m.vz.back();

      // set name, dimensionality, type and connectivity
      ostringstream name;
      name << NAME;
      if (felemstypes.size()>1)
        name << "_t" << felemstypes[i];
      z.n = name.str();
      z.t = felemstypes[i];

      // set new zone element-node connectivity
      z.e2n.resize(felemsnb[i]);
      for (unsigned j=0, k=0; j<e2n.size(); ++j)
        if (e2n[j].t==felemstypes[i])
          z.e2n[k++].n = e2n[j].n;
    }

  }

  f.close();
}

