
#include "GetPot.h"

#include "mlog.h"

namespace {
extern "C" {
#include "ext/ngspice/src/include/ngspice/sharedspice.h"
}
}


using namespace std;


int main(int argc, char **argv)
{
  // user options
  GetPot o(argc,argv);
#if 0
  if (o.search(2,"--help","-h") || o.size()<2) {
    cout << "\n"
         "Adjust working/counter electrodes (WE/CE) metal potentials.\n"
         "(* assumes zero Ohmic drop)\n"
         "Usage (" << o[0] << "):\n"
         "  --help,-h:        this help\n"
         "  --case,-c [str]:  MITReM testcase file:label\n"
         "  -v [str]:         display every iteration (default: no)\n"
         "  -o [str]:         log file to write to (default: none)\n"
         "  -@:               current-fit complementary electrode (default: no)\n"
         "\n"
         "  -Niter [int]:     Newton method max. iterations (default: 1000)\n"
         "  -NdVL2 [real]:    ... max. error L2(dV) norm (default: 1.e-24)\n"
         "  -NdVmax [real]:   ... max. dV update (default: 1.)\n"
         "  -eps [real]:      numerical derivative perturbation (default: 1.e-8)\n"
         "\n"
         "  -I[we|ce] [real[:real:real]]:  WE/CE *intended current (default 0.)\n"
         "  -J[we|ce] [real[:real:real]]:  ... *intended current density (default 0.)\n"
         "  -V[we|ce] [real[:real:real]]:  ... *metal potential value (default 0.)\n"
         "  -R[we|ce] [str:...]:           ... electrode reactions (default: ':' (all))\n"
         "  -A[we|ce] [real]:              ... areas (default: 1.)\n"
         "  (* ranges [min:step:max] only allowed on one parameter)\n"
         "\n"
         "  Specify these parameters to adjust for:\n"
         "  - metal pot. difference:  -V[ce,we] -R[ce,we] -A[ce,we]\n"
         "  - WE or CE current:       -I[we|ce] -R[we|ce] -A[we|ce]\n"
         "  - WE or CE c. density:    -J[we|ce] -R[we|ce]\n"
         << endl;
    return 0;
  }
#endif



  mlog::facility::channel(mlog::DEBUG,new mlog::LogChannelStream(std::cout));


  struct spiceinterface_t{

    /*
     * sending output from stdout, stderr to caller
     * char* string to be sent to caller output
     * int   identification number of calling ngspice shared lib
     * void* return pointer received from caller, e.g. pointer to object having sent the request
     */
    static int printfcn(char* msg, int id, void*) {
      istringstream is(msg);
      string tag;
      is >> tag;
      tag=="stderr"? mlog::warn() << "spice[" << id << "]: " << is.str().substr(tag.length()) << mlog::nl :
      tag=="stdout"? mlog::info() << "spice[" << id << "]: " << is.str().substr(tag.length()) << mlog::nl :
                     mlog::info() << "spice[" << id << "]: " << msg << mlog::nl;
      return 0;
    }

    /*
     * sending simulation status to caller
     * char* simulation status and value (in percent) to be sent to caller
     * int   identification number of calling ngspice shared lib
     * void* return pointer received from caller
     */
    static int statusfcn(char*, int, void*) {
      mlog::info() << "statusfcn" << mlog::nl;
      return 0;
    }

    /*
     * asking for controlled exit
     * int   exit status
     * bool  if true: immediate unloading dll, if false: just set flag, unload is done when function has returned
     * bool  if true: exit upon 'quit', if false: exit due to ngspice.dll error
     * int   identification number of calling ngspice shared lib
     * void* return pointer received from caller
     */
    static int ngspiceexit(int, bool, bool, int, void*) {
      mlog::info() << "ngspiceexit" << mlog::nl;
      return 0;
    }

    /*
     * send back actual vector data
     * vecvaluesall* pointer to array of structs containing actual values from all vectors
     * int           number of structs (one per vector)
     * int           identification number of calling ngspice shared lib
     * void*         return pointer received from caller
     */
    static int sdata(pvecvaluesall, int, int, void*) {
      mlog::info() << "sdata" << mlog::nl;
      return 0;
    }

    /*
     * send back initialization vector data
     * vecinfoall* pointer to array of structs containing data from all vectors right after initialization
     * int         identification number of calling ngspice shared lib
     * void*       return pointer received from caller
     */
    static int sinitdata(pvecinfoall, int, void*) {
      mlog::info() << "sinitdata" << mlog::nl;
      return 0;
    }

    /*
     * indicate if background thread is running
     * bool        true if background thread is running
     * int         identification number of calling ngspice shared lib
     * void*       return pointer received from caller
     */
    static int bgtrun(bool, int, void*) {
      mlog::info() << "bgtrun" << mlog::nl;
      return 0;
    }

  }
  si;


  ngSpice_Init(
    &si.printfcn,
    &si.statusfcn,
    &si.ngspiceexit,
    &si.sdata,
    &si.sinitdata,
    &si.bgtrun,
    static_cast< void* >(&si) );


#if 0
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


#endif
  return 0;
}
