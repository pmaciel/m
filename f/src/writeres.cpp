
#include <fstream>
#include "common.h"

/* outputs information about the non-linear step */
void writeres()
{
  using std::cout;
  using std::endl;

  // write to file
  std::ofstream f(file_log.c_str(),std::ios_base::app);
  f << iter;
  for (int iv=0; iv<Nsys; ++iv)
    f << '\t' << logL2[iv];
  f << endl;
  f.close();

  // write to screen
  cout << "writeres: iteration report..." << endl
       << "  iteration: " << iter << endl
       << "  method: " << (Jacobian==0? "Picard" : (Jacobian==1? "Approximate Newton" : "Newton" )) << endl
       << "  CFL: " << (dtrelax? CFL : 0.) << endl
       << "  Qin=" << Qin << " Qout=" << Qout << " delQ=" << 100.*(Qin+Qout)/(Qin+epsilon) << " %Qin" << endl
       << "  " << "v"               << '\t' << "L1"      << '\t' << "L2"      << '\t' << "Li"      << endl;
  for (int iv=0; iv<Nsys; ++iv)
  cout << "  " <<  m_vars_label[iv] << '\t' << logL1[iv] << '\t' << logL2[iv] << '\t' << logLi[iv] << endl;
  cout << "writeres: iteration report." << endl;
}

