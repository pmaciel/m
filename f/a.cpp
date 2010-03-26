
#include "ext/GetPot.h"
#include "MITReM.h"
#include "ElementMatrixAssembler.h"


double current(MITReM &M, const double &A, const double &V)
{
  double I=0.;
  for (unsigned r=0; r<M.getNElecReactions(); ++r)
    I += A*F_CONST*M.getElecReactionNElectrons(r)*M.calcElecReactionRate(r,V);
  return I;
}


int main(int argc, char **argv)
{
  using namespace std;

  // user options
  GetPot o(argc,argv);
  if (o.search(2,"--help", "-h") || o.size()<2) {
    cout << endl
         << "Adjuster for electrode metal potentials. Assumes:" << endl
         << "  * all electrode reactions are active" << endl
         << "  * zero Ohmic drop" << endl
         << endl
         << "Usage (" << o[0] << "):" << endl
         << "  --help,-h:         this help" << endl
         << "  --case,-c [str]:   MITReM testcase file:label"   << endl
         << endl
         << "  Newton method:" << endl
         << "  -Niter [int]:      Newton method maximum number of iterations (default: 100)" << endl
         << "  -Numax [real]:     ... solution update maximum (default: 1.)" << endl
         << endl
         << "  Adjust for current (either working or counter electrode):" << endl
         << "  -Iwe,-Ice [real]:  intended current (default 0.)" << endl
         << "  -Vwe,-Vce [real]:  m. potential starting condition (default 0.)" << endl
         << "  -Awe,-Ace [real]:  areas (default: 1.)" << endl
         << endl
         << "  Adjust for metal potential difference:" << endl
         << "  -Vwe,-Vce [real]:  intended m. potential difference (and starting condition)" << endl
         << "  -Awe,-Ace [real]:  areas (default: 1.)" << endl
         << endl;
    return 0;
  }
  const string eccase = o.follow("",2,"--case","-c");
  const unsigned Niter = o.follow(100,"-Niter");
  const double   Numax = o.follow(1.,"-Numax");
  const double Vwe = o.follow(0.,"-Vwe"), Vce = o.follow(0.,"-Vce"),
               Iwe = o.follow(0.,"-Iwe"), Ice = o.follow(0.,"-Ice"),
               Awe = o.follow(1.,"-Awe"), Ace = o.follow(1.,"-Ace");
  const double eps = 1e-16;
  cout.precision(12);


  cout << "a: setting MITReM..." << endl;

  // build object
  MITReM m(
    eccase.substr(0,eccase.find(":")),
    eccase.find(":")==string::npos? "*" : eccase.substr(eccase.find(":")+1) );

  // set bulk concentrations and initialize mitrem
  vector< double > bulk(m.getNIons(),0.);
  for (unsigned i=0; i<m.getNIons(); ++i)
    bulk[i] = m.getIonInletConcentration(i);
  m.init(&bulk[0], /*U*/ 0., m.getSolutionTemperature(),m.getSolutionDensity());

  cout << "a: setting MITReM." << endl;


  if (o.search("-Iwe") || o.search("-Ice")) {
    const double Ie = (o.search("-Iwe")? Iwe : (o.search("-Ice")? Ice : 0. )),
                 V  = (o.search("-Iwe")? Vwe : (o.search("-Ice")? Vce : 0. )),
                 A  = (o.search("-Iwe")? Awe : (o.search("-Ice")? Ace : 0. ));
    cout << "a: forced current (Ie [A] = " << Ie << ")..." << endl;

    double dV = 0.;
    for (unsigned i=1; i<=Niter; ++i) {
      const double I =  current(m,A,V+dV),
                  dI = (current(m,A,V+dV+eps)-I)/eps,
                   x = -(I-Ie)/dI;
      dV += (x<0.? -1.:1.) * min(Numax,abs(x));
      cout << "a:  i=" << i << "  I [A] = " << current(m,A,V+dV) << endl;
    }

    cout << "a: adjusted (I-Ie)/V [A/V]: " << current(m,A,V+dV)-Ie << " / " << V+dV << endl;
    cout << "a: forced current." << endl;
  }
  else {
    cout << "a: forced metal potential difference (Vce-Vwe [V] = " << Vce-Vwe << ")..." << endl;
    cout << "a: original Vwe/Vce [V/V]: " << Vwe << " / " << Vce << endl;

    // Newton method (on current balance)
    double dV = 0.;  // shift initial guess
    for (unsigned i=1; i<=Niter; ++i) {

      double  I = current(m,Awe,Vwe+dV)     + current(m,Ace,Vce+dV),
             dI = current(m,Awe,Vwe+dV+eps) + current(m,Ace,Vce+dV+eps);
      dI = (dI-I)/eps;

      dV += -I/dI;
      cout << "a:  i=" << i << "  dV/|ddV| [V/V] = " << dV << " / " << abs(-I/dI) << endl;

    }

    cout << "a: adjusted Vwe/Vce [V/V]: " << Vwe+dV << " / " << Vce+dV << endl;
    cout << "a: forced metal potential difference." << endl;
  }

  return 0;
}

