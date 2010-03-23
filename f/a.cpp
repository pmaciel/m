
#include "ext/GetPot.h"
#include "MITReM.h"
#include "ElementMatrixAssembler.h"

int main(int argc, char **argv)
{
  using namespace std;

  // user options
  GetPot o(argc,argv);
  if (o.search(2,"--help", "-h")) {
    cout << endl
         << "Usage (" << o[0] << "):" << endl
         << "  --help, -h: this help" << endl
         << "  --case, -c: testcase file:label (default: \"mutech.ec.xml\")"   << endl
         << "  -Vwe:       metal potential at the working electrode (default: 0.)" << endl
         << "  -Vwe:       ... counter electrode (default: 1.)"                    << endl
         << "  -Awe:       area of the working electrode (default: 1.)"            << endl
         << "  -Awe:       ... counter electrode (default: 1.)"                    << endl
         << "  -Niter:     maximum number of iterations to converge for (default: 100)" << endl
         << endl;
    return 0;
  }
  const string eccase = o.follow("mutech.ec.xml",2,"--case","-c");
  double Vwe = o.follow(0.,"-Vwe"),
         Vce = o.follow(1.,"-Vce");
  const double Awe = o.follow(1.,"-Awe"),
               Ace = o.follow(1.,"-Ace");
  const unsigned Niter = o.follow(100,"-Niter");


  cout << "a: setting MITReM..." << endl;
  MITReM m(
    eccase.substr(0,eccase.find(":")),
    eccase.find(":")==string::npos? "*" : eccase.substr(eccase.find(":")+1) );
  cout << "a: setting MITReM." << endl;


  cout << "a: solve..." << endl;
  m.correctVForPotentialDifference(Vwe,Vce,Awe,Ace,Niter);
  cout << "a: solve." << endl;


  return 0;
}

