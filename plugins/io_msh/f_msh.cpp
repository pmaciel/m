
#include "mfactory.h"
#include "f_msh.h"
#include <sstream>

using namespace std;
using namespace m;


Register< mfinput,f_msh > mf_msh(".msh","Gmsh mesh input format");


namespace aux {


enum        msh_section_t     { NO_SECTION=0,  MESHFORMAT,   NODES,   ELEMENTS,   PERIODIC,   PHYSICALNAMES,   NODEDATA,   ELEMENTDATA,   INTERPOLATIONSCHEME,   ALL_SECTIONS };
const char* msh_section_n[] = { "",           "MeshFormat", "Nodes", "Elements", "Periodic", "PhysicalNames", "NodeData", "ElementData", "InterpolationScheme", "" };


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


template< int S=NO_SECTION >
struct msh_section_reader {
  msh_section_reader() : sbegin(msh_section_n[S]), send("End"+sbegin) {}
  const std::string
    sbegin,
    send;
};


template< int E=0 >
struct msh_element_reader {
  msh_element_reader() : nnodes(element_n[E]) {}
  const unsigned nnodes;
};


}


void f_msh::read(GetPot& o,mmesh& m)
{
  const string filename(o.get(o.inc_cursor(),""));
  ifstream f(filename.c_str());
  if (!f) {
    cerr << "error accessing file: \"" << filename << "\"" << endl;
    throw 42;
  }


  // here


  f.close();
}

