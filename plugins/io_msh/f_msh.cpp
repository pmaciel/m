#include <fstream>
#include <sstream>

#include <boost/lexical_cast.hpp>

#include "mfactory.h"
#include "f_msh.h"

using namespace std;
using namespace m;


Register< mfinput,f_msh > mf_msh(".msh","Gmsh mesh input format");


namespace aux {


enum        msh_section_t     { NO_SECTION=0,  MESHFORMAT,   NODES,   ELEMENTS,   PERIODIC,   PHYSICALNAMES,   NODEDATA,   ELEMENTDATA,   INTERPOLATIONSCHEME,   COMMENTS,   ALL_SECTIONS };
const char* msh_section_n[] = { "",           "MeshFormat", "Nodes", "Elements", "Periodic", "PhysicalNames", "NodeData", "ElementData", "InterpolationScheme", "Comments", "" };


const int element_t[] = {
  0,
  FELINESEG,        //  1:   2-node line
  FETRIANGLE,       //  2:   3-node triangle
  FEQUADRILATERAL,  //  3:   4-node quadrangle
  FETETRAHEDRON,    //  4:   4-node tetrahedron
  FEBRICK,          //  5:   8-node hexahedron
  PRISM3,           //  6:   6-node prism
  PYRAMID4,         //  7:   5-node pyramid
  0,                //  8:   3-node second order line (2 nodes associated with the vertices and 1 with the edge)
  0,                //  9:   6-node second order triangle (3 nodes associated with the vertices and 3 with the edges)
  0,                // 10:   9-node second order quadrangle (4 nodes associated with the vertices, 4 with the edges and 1 with the face)
  0,                // 11:  10-node second order tetrahedron (4 nodes associated with the vertices and 6 with the edges)
  0,                // 12:  27-node second order hexahedron (8 nodes associated with the vertices, 12 with the edges, 6 with the faces and 1 with the volume)
  0,                // 13:  18-node second order prism (6 nodes associated with the vertices, 9 with the edges and 3 with the quadrangular faces)
  0,                // 14:  14-node second order pyramid (5 nodes associated with the vertices, 8 with the edges and 1 with the quadrangular face)
  ORDERED,          // 15:   1-node point
  0,                // 16:   8-node second order quadrangle (4 nodes associated with the vertices and 4 with the edges)
  0,                // 17:  20-node second order hexahedron (8 nodes associated with the vertices and 12 with the edges)
  0,                // 18:  15-node second order prism (6 nodes associated with the vertices and 9 with the edges)
  0,                // 19:  13-node second order pyramid (5 nodes associated with the vertices and 8 with the edges)
  0,                // 20:   9-node third order incomplete triangle (3 nodes associated with the vertices, 6 with the edges)
  0,                // 21:  10-node third order triangle (3 nodes associated with the vertices, 6 with the edges, 1 with the face)
  0,                // 22:  12-node fourth order incomplete triangle (3 nodes associated with the vertices, 9 with the edges)
  0,                // 23:  15-node fourth order triangle (3 nodes associated with the vertices, 9 with the edges, 3 with the face)
  0,                // 24:  15-node fifth order incomplete triangle (3 nodes associated with the vertices, 12 with the edges)
  0,                // 25:  21-node fifth order complete triangle (3 nodes associated with the vertices, 12 with the edges, 6 with the face)
  0,                // 26:   4-node third order edge (2 nodes associated with the vertices, 2 internal to the edge)
  0,                // 27:   5-node fourth order edge (2 nodes associated with the vertices, 3 internal to the edge)
  0,                // 28:   6-node fifth order edge (2 nodes associated with the vertices, 4 internal to the edge)
  0,                // 29:  20-node third order tetrahedron (4 nodes associated with the vertices, 12 with the edges, 4 with the faces)
  0,                // 30:  35-node fourth order tetrahedron (4 nodes associated with the vertices, 18 with the edges, 12 with the faces, 1 in the volume)
  0,                // 31:  56-node fifth order tetrahedron (4 nodes associated with the vertices, 24 with the edges, 24 with the faces, 4 in the volume)
      0,0,0,0,0,0,0,0,  // (not assigned)
  0,0,0,0,0,0,0,0,0,0,
  0,0,0,0,0,0,0,0,0,0,
  0,0,0,0,0,0,0,0,0,0,
  0,0,0,0,0,0,0,0,0,0,
  0,0,0,0,0,0,0,0,0,0,
  0,0,
  0,                // 92:  64-node third order hexahedron (8 nodes associated with the vertices, 24 with the edges, 24 with the faces, 8 in the volume)
  0                 // 93: 125-node fourth order hexahedron (8 nodes associated with the vertices, 36 with the edges, 54 with the faces, 27 in the volume)
};


const unsigned element_n[] = {
    0,  2,  3,  4,  4,  8,  6,  5,  3,  6,
    9, 10, 27, 18, 14,  1,  8, 20, 15, 13,
    9, 10, 12, 15, 15, 21,  4,  5,  6, 20,
   35, 56,
        0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,0,0,
    0,0,
   64, 125
};


template< int E=0 >
struct msh_element_reader {
  msh_element_reader() : nnodes(element_n[E]) {}
  const unsigned nnodes;
};


}  // namespace aux


void f_msh::read(GetPot& o,mmesh& m)
{
  const string filename(o.get(o.inc_cursor(),""));
  ifstream f(filename.c_str());
  if (!f) {
    cerr << "error accessing file: \"" << filename << "\"" << endl;
    throw 42;
  }


  // read per section (magic!)
  m.vn.clear();
  m.vv.clear();
  m.vz.clear();
  while (f) {
    string s;
    f >> s;
    const string end("$End"+(s.length()? s.substr(1):s));
    s=="$MeshFormat"?          readMeshFormat          (end,f,m) :
    s=="$Nodes"?               readNodes               (end,f,m) :
    s=="$Elements"?            readElements            (end,f,m) :
    s=="$Periodic"?            readPeriodic            (end,f,m) :
    s=="$PhysicalNames"?       readPhysicalNames       (end,f,m) :
    s=="$NodeData"?            readNodeData            (end,f,m) :
    s=="$ElementData"?         readElementData         (end,f,m) :
    s=="$InterpolationScheme"? readInterpolationScheme (end,f,m) :
    s=="$Comments"?            readComments(end,f) :
                               readError(f);
  }


  f.close();
}


void f_msh::readMeshFormat(const string& end, ifstream& f, mmesh& m)
{
  double versionnumber(-1.);
  int filetype(-1), datasize(-1);
  string stop;
  f? (f>>versionnumber)? (f>>filetype)? (f>>datasize)? (f>>stop)? : f:f:f:f:f;
  if ( ((versionnumber!=2.2) && (cerr << "msh: warning: version-number!=2.2"  << endl)) ||
       ((filetype!=0)        && (cerr << "msh: warning: file-type!=0"         << endl)) ||
       ((datasize!=8)        && (cerr << "msh: warning: data-size!=8"         << endl)) ||
       ((stop!=end)          && (cerr << "msh: warning: section not finished" << endl)) )
    f.close();
}


void f_msh::readNodes(const string& end, ifstream& f, mmesh& m)
{
  unsigned numberofnodes(0), n(1);
  f? (f>>numberofnodes) : f;

  m.vn = vector< string >(3,"x");  m.vn[1]="y";  m.vn[2]="z";
  m.vv.assign(3,vector< double >(numberofnodes,0.));

  while (f && n>=1 && n<=numberofnodes) {
    f >> n;
    f >> m.vv[0][n-1] >> m.vv[1][n-1] >> m.vv[2][n-1];
  }
}


void f_msh::readElements(const string& end, ifstream& f, mmesh& m)
{
  unsigned numberofelements(0), elmnumber(1), elmtype(0);
  f? (f>>numberofelements) : f;

  //FIXME does not support mixed element types
  m.vz.push_back(m::mzone());
  mzone& z = m.vz.back();
  z.n = "Zone"+boost::lexical_cast< string >(m.vz.size());
  for (; (f>>elmnumber>>elmtype) && (elmnumber>=1) && (elmnumber<=numberofelements); ) {

    z.t = static_cast< mtype >(aux::element_t[elmtype]);
    if (numberofelements && !(z.e2n[0].n.size())) {
      melem e;
      e.n.assign(aux::element_n[elmtype],0);
      z.e2n.assign(numberofelements,e);
    }

    vector< unsigned > &n = z.e2n[elmnumber-1].n;
    n.assign(aux::element_n[elmtype],0);
    for (size_t i=0; i<n.size(); ++i)
      f >> n[i];

  }
}


void f_msh::readPeriodic(const string& end, ifstream& f, const mmesh& m)
{
}


void f_msh::readPhysicalNames(const string& end, ifstream& f, mmesh& m)
{
}


void f_msh::readNodeData(const string& end, ifstream& f, mmesh& m)
{
}


void f_msh::readElementData(const string& end, ifstream& f, const mmesh& m)
{
}


void f_msh::readInterpolationScheme(const string& end, ifstream& f, const mmesh& m)
{
}


void f_msh::readComments(const string& end, ifstream& f)
{
  string s;
  while ((f>>s) && (s!=end)) {};
}


void f_msh::readError(ifstream& f)
{
  f.close();
}
