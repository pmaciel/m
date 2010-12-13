
#include "ext/GetPot.h"
#include "MITReM.h"
#include "ElementMatrixAssembler.h"


// forward declarations
double cdensity(const MITReM &M, const std::vector< unsigned >& R, const double &V);
std::vector< unsigned > getElecReactions(const MITReM &M, const std::string& s);
std::vector< double > getRangedParam(const std::string& s, double mult=1.);


int main(int argc, char **argv)
{
  using std::cout;
  using std::endl;
  using std::vector;

  // user options
  GetPot o(argc,argv);
  if (o.search(2,"--help","-h") || o.size()<2) {
    cout << endl
         << "Adjust working/counter electrodes (WE/CE) metal potentials." << endl
         << "(* assumes zero Ohmic drop)" << endl
         << "Usage (" << o[0] << "):" << endl
         << "  --help,-h:        this help" << endl
         << "  --case,-c [str]:  MITReM testcase file:label"   << endl
         << "  -v [str]:         display every iteration (default: no)"   << endl
         << "  -o [str]:         log file to write to (default: none)"   << endl
         << "  -@:               current-fit complementary electrode (default: no)" << endl
         << endl
         << "  -Niter [int]:     Newton method max. iterations (default: 1000)" << endl
         << "  -NdVL2 [real]:    ... max. error L2(dV) norm (default: 1.e-24)" << endl
         << "  -NdVmax [real]:   ... max. dV update (default: 1.)" << endl
         << "  -eps [real]:      numerical derivative perturbation (default: 1.e-8)" << endl
         << endl
         << "  -I[we|ce] [real[:real:real]]:  WE/CE *intended current (default 0.)" << endl
         << "  -J[we|ce] [real[:real:real]]:  ... *intended current density (default 0.)" << endl
         << "  -V[we|ce] [real[:real:real]]:  ... *metal potential value (default 0.)" << endl
         << "  -R[we|ce] [str:...]:           ... electrode reactions (default: ':' (all))" << endl
         << "  -A[we|ce] [real]:              ... areas (default: 1.)" << endl
         << "  (* ranges [min:step:max] only allowed on one parameter)" << endl
         << endl
         << "  Specify these parameters to adjust for:" << endl
         << "  - metal pot. difference:  -V[ce,we] -R[ce,we] -A[ce,we]" << endl
         << "  - WE or CE current:       -I[we|ce] -R[we|ce] -A[we|ce]" << endl
         << "  - WE or CE c. density:    -J[we|ce] -R[we|ce]" << endl
         << endl;
    return 0;
  }


  cout << "a: setting MITReM..." << endl;

  // build object
  const std::string eccase = o.follow("",2,"--case","-c");
  MITReM m(
    eccase.substr(0,eccase.find(":")),
    eccase.find(":")==std::string::npos? "*" : eccase.substr(eccase.find(":")+1) );

  // set bulk concentrations and initialize mitrem
  vector< double > bulk(m.getNIons(),0.);
  for (unsigned i=0; i<m.getNIons(); ++i)
    bulk[i] = m.getIonInletConcentration(i);
  m.init(&bulk[0], /*U*/ 0., m.getSolutionTemperature(),m.getSolutionDensity());

  // set electrode reactions
  vector< unsigned >
    Rwe(getElecReactions(m,o.follow(":","-Rwe"))),
    Rce(getElecReactions(m,o.follow(":","-Rce")));

  cout << "a: setting MITReM." << endl;


  // Newton method and other parameters
  const unsigned Niter = o.follow(1000,"-Niter");
  const double dVL2  = o.follow(1.e-24,"-NdVL2"),
               dVmax = o.follow(1.,"-NdVmax"),
               eps = o.follow(1.e-8,"-eps"),
               Awe(o.follow(1.,"-Awe")),
               Ace(o.follow(1.,"-Ace"));
  const bool verbose(o.search("-v"));


  // streams to use
  std::ofstream fout;
  if (o.search("-o")) {
    fout.open(o.follow("a.log","-o").c_str());
    std::ostringstream s;
    s << o[0];
    for (unsigned i=1; i<o.size(); ++i)
      s << ' ' << o[i];
    fout << "command: " << s.str() << endl;
  }
  cout.precision(12);


  // get parameter ranges (avoid storing I[we|ce])
  const vector< double >
    rangeVwe(getRangedParam(o.follow("0.","-Vwe"))),
    rangeVce(getRangedParam(o.follow("0.","-Vce"))),
    rangeJwe(o.search("-Iwe")? getRangedParam(o.follow("0.","-Iwe"),1./Awe)
                             : getRangedParam(o.follow("0.","-Jwe")) ),
    rangeJce(o.search("-Ice")? getRangedParam(o.follow("0.","-Ice"),1./Ace)
                             : getRangedParam(o.follow("0.","-Jce")) );


  // find ranged parameter and step on it
  enum T_RANGEDPARAM { NONE, JWE, JCE, VWE, VCE };
  const T_RANGEDPARAM rangedparam(
    (rangeJwe.size()>1? JWE :
    (rangeJce.size()>1? JCE :
    (rangeVwe.size()>1? VWE :
    (rangeVce.size()>1? VCE : NONE )))) );
  const vector< double >& range(
    (rangedparam==JWE? rangeJwe :
    (rangedparam==JCE? rangeJce :
    (rangedparam==VWE? rangeVwe :
                       rangeVce ))) );
  for (vector< double >::const_iterator p=range.begin(); p!=range.end(); ++p) {


    // set WE/CE parameters
    double
      Jwe(rangedparam==JWE? *p:rangeJwe[0]),
      Jce(rangedparam==JCE? *p:rangeJce[0]),
      Vwe(rangedparam==VWE? *p:rangeVwe[0]),
      Vce(rangedparam==VCE? *p:rangeVce[0]);


    // select method
    if (o.search(4,"-Iwe","-Jwe","-Ice","-Jce")) {
      /*
       * forced current density
       * Newton method to current density value, f(V)=J
       */

      const vector< unsigned >& R = (o.search(2,"-Iwe","-Jwe")? Rwe:Rce);
      const double A(o.search(2,"-Iwe","-Jwe")? Awe:Ace),
                   J(o.search(2,"-Iwe","-Jwe")? Jwe:Jce);
      double       V(o.search(2,"-Iwe","-Jwe")? Vwe:Vce);
      double x(1.);  // (initial) solution update

      for (unsigned i=1; i<=Niter &&
           x*x>dVL2 && !(x!=x || x>1.e100); ++i) {
        const double f  = cdensity(m,R,V),
                     fp = cdensity(m,R,V+eps),
                     df = (fp-f)/eps;
        x = -(f-J)/df;
        V += (x<0.? -1.:1.) * std::min(dVmax,std::abs(x));
        if (verbose)
          cout << "a:  i=" << i << "  L2=" << x*x
               << "  J[A.m-2]=" <<   cdensity(m,R,V)
               <<     "  I[A]=" << A*cdensity(m,R,V) << endl;
      }
      if (o.search(2,"-Iwe","-Jwe")) { Vwe=V; Vce=1./0.; }
      else                           { Vce=V; Vwe=1./0.; }

      if (o.search("-@")) {
        /*
         * fit complementary electrode
         */

        const vector< unsigned >& _R = (o.search(2,"-Iwe","-Jwe")? Rce:Rwe);
        const double _A(o.search(2,"-Iwe","-Jwe")? Ace:Awe),
                     _J(-J*A/_A);  // necessary J for a balanced I
        double       _V(-V);       // complementary V starts symmetric
        x = 1.;  // (initial) solution update

        for (unsigned i=1; i<=Niter &&
             x*x>dVL2 && !(x!=x || x>1.e100); ++i) {
          const double f  = cdensity(m,_R,_V),
                       fp = cdensity(m,_R,_V+eps),
                       df = (fp-f)/eps;
          x = -(f-_J)/df;
          _V += (x<0.? -1.:1.) * std::min(dVmax,std::abs(x));
          if (verbose)
            cout << "a:  i=" << i << "  L2=" << x*x
                 << "  J[A.m-2]=" <<    cdensity(m,_R,_V)
                 <<     "  I[A]=" << _A*cdensity(m,_R,_V) << endl;
        }
        if (o.search(2,"-Iwe","-Jwe")) { Vce=_V; }
        else                           { Vwe=_V; }
      }

    }
    else {
      /*
       * forced metal potential difference
       * Newton method to current balance, f(V)=sum(I)=0.
       */

      double x(1.);  // solution update
      for (unsigned i=1; i<=Niter &&
           x*x>dVL2 && !(x!=x || x>1.e100); ++i) {
        const double f  = Awe*cdensity(m,Rwe,Vwe)     + Ace*cdensity(m,Rce,Vce),
                     fp = Awe*cdensity(m,Rwe,Vwe+eps) + Ace*cdensity(m,Rce,Vce+eps),
                     df = (fp-f)/eps;
        x = -(f-0.)/df;
        Vwe += (x<0.? -1.:1.) * std::min(dVmax,std::abs(x));
        Vce += (x<0.? -1.:1.) * std::min(dVmax,std::abs(x));
        if (verbose)
          cout << "a:  i=" << i << "  L2=" << x*x << endl;
      }

    }


    // display/write adjustment message
    {
      const double _Jwe(cdensity(m,Rwe,Vwe)),
                   _Jce(cdensity(m,Rce,Vce));
      {
        cout << "a: adjustment: "
             << (o.search(2,"-Iwe","-Jwe")? "(Jwe[A.m-2]=" :
                (o.search(2,"-Ice","-Jce")? "(Jce[A.m-2]=" : "(Vce-Vwe[V]=" ))
             << (o.search(2,"-Iwe","-Jwe")? Jwe :
                (o.search(2,"-Ice","-Jce")? Jce : Vce-Vwe )) << ')'
             << "  Vce:"  << Vce << "  Jce:"  << _Jce << "  Ice:"  << Ace*_Jce
             << "  Vwe:"  << Vwe << "  Jwe:"  << _Jwe << "  Iwe:"  << Awe*_Jwe
             << "  (Vce-Vwe):"  << (Vce-Vwe)
             << "  |Ice+Iwe|:"  << std::abs(Ace*_Jce + Awe*_Jwe)
             << endl;
      }
      if (fout) {
        fout << (o.search(2,"-Iwe","-Jwe")? "(Jwe[A.m-2]=" :
                (o.search(2,"-Ice","-Jce")? "(Jce[A.m-2]=" : "(Vce-Vwe[V]=" ))
             << (o.search(2,"-Iwe","-Jwe")? Jwe :
                (o.search(2,"-Ice","-Jce")? Jce : Vce-Vwe )) << ')'
             << "\tVce\t" << Vce << "\tJce\t" << _Jce << "\tIce\t" << Ace*_Jce
             << "\tVwe\t" << Vwe << "\tJwe\t" << _Jwe << "\tIwe\t" << Awe*_Jwe
             << "\t(Vce-Vwe)\t" << (Vce-Vwe)
             << "\t|Ice+Iwe|\t" << std::abs(Ace*_Jce + Awe*_Jwe)
             << endl;
      }
    }


  }  // step on ranged parameter


  return 0;
}


double cdensity(const MITReM &M, const std::vector< unsigned >& R, const double &V)
{
  double J(0.);
  for (std::vector< unsigned >::const_iterator r=R.begin(); r!=R.end(); ++r)
    J += F_CONST*M.getElecReactionNElectrons(*r)*M.calcElecReactionRate(*r,V);
  return J;
}


template< typename T >
std::vector< T > getVValues(const std::string& l)
{
  using std::string;
  using std::vector;
  using std::istringstream;

  // split string by ':'
  vector< T > r;
  if (l.length() && l!=":") {
    string::size_type p1 = 0,
                      p2 = 0;
    while (p2!=string::npos) {
      p2 = l.find(":",p1);
      istringstream ss( l.substr(p1,(p2==string::npos? p2:p2-p1)) );
      r.push_back(T());
      ss >> r.back();
      p1 = p2+1;
    }
  }
  return r;
}


std::vector< unsigned > getElecReactions(const MITReM &M, const std::string& s)
{
  using std::vector;
  using std::string;
  vector< unsigned > r;

  const vector< string > vs(getVValues< std::string >(s));
  if (!vs.size()) {
    // initialize reactions if all are to be used
    r.assign(M.getNElecReactions(),0);
    for (unsigned i=1; i<r.size(); ++i)
      r[i] = r[i-1]+1;
  }
  else {
    // otherwise search reaction labels to build vector of indices
    for (vector< string >::const_iterator l=vs.begin(); l!=vs.end(); ++l)
      for (unsigned i=0; i<M.getNElecReactions(); ++i)
        if (*l==M.getElecReactionLabel(i)) {
          r.push_back(i);
          break;
        }
    if (vs.size()!=r.size() && vs.size()) {
      std::cerr << "a: some reactions were not found!" << std::endl;
      throw 42;
    }
  }
  return r;
}


std::vector< double > getRangedParam(const std::string& s, double mult)
{
  using std::vector;
  vector< double > param;

  const vector< double > vd(getVValues< double >(s));
  if      (vd.size()==1) { param.assign(1,vd[0]); }
  else if (vd.size()==3) {
    const double step((vd[1]/vd[1]!=vd[1]/vd[1]) || (vd[1]/vd[1]>1.e100)? 0. : vd[1]),
                 _min(step>0.? std::min(vd[0],vd[2]):std::max(vd[0],vd[2])),
                 _max(step<0.? std::min(vd[0],vd[2]):std::max(vd[0],vd[2]));
    if (std::fabs(step)<1.e-12) {
      std::cerr << "a: invalid range step (|step|<1.e-12)!" << std::endl;
      throw 42;
    }
    param.assign(1,_min);
    while (step>0.? param.back()<_max : param.back()>_max)
      param.push_back(param.back()+step);
    if (std::fabs(param.back())>std::fabs(_max))
      param.back() = _max;
  }
  else {
    std::cerr << "a: incorrect range specified, 1 or 3 parameters expected!" << std::endl;
    throw 42;
  }
  for (vector< double >::iterator p=param.begin(); p!=param.end(); ++p)
    *p *= mult;

  return param;
}
