
#include "mfactory.h"
#include "f_cgns.h"
#include <sstream>

using namespace std;
using namespace m;


Register< mfinput, f_cgns > mf_cgns1(".cgns","CGNS input format (MARIN, GridPro and Hexpress variants)");
Register< mfoutput,f_cgns > mf_cgns2(".cgns","CGNS output format (MARIN variant)");


void f_cgns::read(GetPot& o,mmesh& m)
{
  const string filename(o.get(o.inc_cursor(),""));
  ifstream f(filename.c_str());
  if (!f) {
    cerr << "error accessing file: \"" << filename << "\"" << endl;
    throw 42;
  }


  // ideas accepted


  f.close();
}


void f_cgns::write(GetPot& o, const mmesh& m)
{
  const string fn(o.get(o.inc_cursor(),""));
  ofstream f(fn.c_str());
  if (!f) {
    cerr << "error accessing file: \"" << fn << "\"" << endl;
    throw 42;
  }


  // (even more) ideas accepted


  f.close();
}

