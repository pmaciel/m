
#include "mfactory.h"
//include "f_plt.h"
#include "t_solution.h"

using namespace std;
using namespace m;


Register< mtransform,t_solution > mt_solution("-tsol","[str]: Tecplot (nodal) solution file");


namespace aux {
  string extension(const string& fn)
  {
    const string::size_type idx = fn.rfind('.');
    return fn.substr(idx!=string::npos? idx:0);
  }
}


void t_solution::transform(GetPot& o, mmesh& mold)
{
  const string fn = o.get(o.inc_cursor(),"");
  mmesh mnew;


  cout << "read: \"" << fn << "\"..." << endl;
  {
    const int argc = 2;
    char*     argv[] = { (char*) "", const_cast< char* >(fn.c_str()) };
    GetPot o2(argc,argv);
    const string key(aux::extension(fn));
    mfinput* p = mfactory< mfinput >::instance()->Create(key);
    p->read(o2,mnew);
    delete p;
  }
  cout << "read: \"" << fn << "\"." << endl;


  // some checks
  if (mnew.d()!=mold.d())
    cout << "warn: different dimensions, elements might become damaged!" << endl;
  if (mnew.n()!=mold.n()) {
    cout << "warn: different number of nodes, not applying!" << endl;
    return;
  }
  if (!mnew.v()) {
    cout << "warn: no variables found, not applying!" << endl;
    return;
  }


  cout << "info: replace point cloud..." << endl;
  mold.vn.swap(mnew.vn);
  mold.vv.swap(mnew.vv);
  cout << "info: replace point cloud." << endl;
}


