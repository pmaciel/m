
#include "mfactory.h"
#include "t_solution.h"

using namespace m;


Register< mtransform,t_solution > mt_solution(2,"-tsol","[str] read solution from another file (substitution)",
                                                "-tapp","[str] read solution from another file (append)");


void t_solution::transform(GetPot& o, mmesh& mold, const XMLNode& x)
{
  using namespace std;
  const string k  = o[o.get_cursor()];
  const string fn = o.get(o.inc_cursor(),"");
  mmesh mnew;


  cout << "read: \"" << fn << "\"..." << endl;
  {
    //FXIME dirty hack creating new command line options, to be removed
    const int argc = 2;
    char*     argv[] = { (char*) "", const_cast< char* >(fn.c_str()) };
    GetPot o2(argc,argv);

    mfinput* p = mfactory< mfinput >::instance()->Create(utils::get_file_extension(fn));
    p->read(o2,mnew,x);
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


  if (k=="-tsol") {
    cout << "info: replace point cloud..." << endl;
    mold.vn.swap(mnew.vn);
    mold.vv.swap(mnew.vv);
    cout << "info: replace point cloud." << endl;
  }
  else if (k=="-tapp") {
    cout << "info: append point cloud..." << endl;
    // append variables with different names
    for (unsigned i=0; i<mnew.v(); ++i)
      if (!count(mold.vn.begin(),mold.vn.end(),mnew.vn[i])) {
        cout << "info: appending variable: \"" << mnew.vn[i] << '"' << endl;
        mold.vn.push_back(mnew.vn[i]);
        mold.vv.push_back(mnew.vv[i]);
      }
    cout << "info: append point cloud." << endl;
  }
}

