
#include "mfactory.h"
#include "t_info.h"

using namespace std;
using namespace m;


Register< mtransform,t_info > mt_info("-ti","current mesh information");


void t_info::transform(GetPot& o, mmesh& m)
{
  using std::cout;
  using std::endl;

  cout << "  mesh:"
       << "  d=" << m.d()
       << "  n=" << m.n()
       << "  v=" << m.v()
       << "  z=" << m.z() << endl;
  for (unsigned t=0; t<m.v(); ++t)
    cout << "  variable[" << t << "]: \"" << m.vn[t] << "\"" << endl;
  for (unsigned t=0; t<m.z(); ++t)
    cout << "  zone[" << t << "]:"
         << " d=" << m.d(t)
         << " e=" << m.e(t)
         << " n=\"" << m.vz[t].n << '"'
         << " t=" << (m.vz[t].t==ORDERED?         "ORDERED"         :
                     (m.vz[t].t==FELINESEG?       "FELINESEG"       :
                     (m.vz[t].t==FETRIANGLE?      "FETRIANGLE"      :
                     (m.vz[t].t==FEQUADRILATERAL? "FEQUADRILATERAL" :
                     (m.vz[t].t==FETETRAHEDRON?   "FETETRAHEDRON"   :
                     (m.vz[t].t==FEBRICK?         "FEBRICK"         :
                     (m.vz[t].t==FEPOLYGON?       "FEPOLYGON"       :
                     (m.vz[t].t==FEPOLYHEDRON?    "FEPOLYHEDRON"    :
                     (m.vz[t].t==PRISM3?          "PRISM3"          :
                     (m.vz[t].t==PYRAMID4?        "PYRAMID4"        :
                                                  "(unknown)" ))))))))))

         << endl;
}

