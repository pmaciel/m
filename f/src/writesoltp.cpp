
/* writes solution in tecplot format */

#include <memory>
#include "common.h"

void writesoltp(const std::string& outfile)
{
  std::ostringstream fn;
  fn << outfile << '.' << iter << ".plt";
  std::cout << "writesoltp: writing \"" << fn.str() << "\"..." << std::endl;

  // copy variables: pressure (iv=0) is multiplied by density (rho)
  for (int n=0; n<Nnode; ++n)
    M.vv[Ndim+0][n] = No_W[0][n]*rho;
  for (int iv=1; iv<Nsys; ++iv)
    for (int n=0; n<Nnode; ++n)
      M.vv[Ndim+iv][n] = No_W[iv][n];

  // if wall distance is used, append it as the last variable (remove later)
  if (walldist) {
    M.vn.push_back("wd");
    M.vv.resize(M.vv.size()+1);
    M.vv.back().swap(No_wd); // no copy, better
  }

  // write the file
  e2n.swap(M.vz[0].e2n);  // (hack)
  {
    std::auto_ptr< m::mfoutput > p(m::mfactory< m::mfoutput  >::instance()->Create(".plt"));
    const std::string f = fn.str();
    char* argv[] = { (char*) "", const_cast< char* >(f.c_str()) };
    GetPot o2(2,argv);
    p->write(o2,M);
  }
  e2n.swap(M.vz[0].e2n);  // (hack back)

  // if wall distance is used, remove it
  if (walldist) {
    M.vn.pop_back();
    M.vv.back().swap(No_wd); // no copy, better
    M.vv.pop_back();
  }

  std::cout << "writesoltp: writing." << std::endl;
}
