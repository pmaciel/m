#ifndef mfunction_h
#define mfunction_h

#include <cmath>
#include <cassert>
#include <vector>
#include <string>
#include "ext/fparser.hh"


namespace m {


// tool for evaluating one user-evaluatable function
class mfunction {

 public:
  // constructor with function and variables definitions
  mfunction(const std::string & _fs, const std::vector< std::string >& _vvs)
  {
    // add constants
    f.AddConstant("pi",M_PI);

    using namespace std;
    assert(_vvs.size());
    string vars(_vvs[0]);
    for (int i=1; i<(int) _vvs.size(); ++i)
      vars += ',' + _vvs[i];
    clog << "m::mfunction::f(" << vars << ") = \"" << _fs << '"' << endl;

    const int r = f.Parse(_fs,vars);
    if (r>=0) {
      cerr << "m::mfunction::error:" << string(r+1+vars.length(),' ') << '^' << endl
           << f.ErrorMsg() << endl;
      throw 42;
    }
  }

  // evaluate (up to Nr-th) functions at given variables vector address
  double eval(double* v) {
    return f.Eval(v);
  }

 private:
  // function
  FunctionParser f;

};


// tool for evaluating a vector of user-evaluatable functions
class mvfunction {

 public:
  // constructor with functions and variables definitions
  mvfunction(const std::vector< std::string >& _vfs, const std::vector< std::string >& _vvs) :
    vf(_vfs.size())
  {
    using namespace std;
    assert(_vvs.size());
    string vars(_vvs[0]);
    for (int i=1; i<(int) _vvs.size(); ++i)
      vars += ',' + _vvs[i];
    clog << "m::mvfunction::vars: \"" << vars << '"' << endl;

    assert(_vfs.size());
    for (int i=0; i<(int) _vfs.size(); ++i) {
      const int r = vf[i].Parse(_vfs[i],vars);
      if (r>=0) {
        cerr << string(r+7, ' ') << '^' << endl
             << vf[i].ErrorMsg() << endl;
        throw 42;
      }
      clog << "m::mvfunction::funs: \"" << _vfs[i] << '"' << endl;
    }
  }

  // evaluate (up to Nr-th) functions at given variables vector address
  void eval(int Nr, double* r, double* v) {
    assert(Nr<=(int) vf.size());
    for (int i=0; i<Nr; ++i)
      r[i] = vf[i].Eval(v);
  }

 private:
  // vector of functions
  std::vector< FunctionParser > vf;

};


}  // namespace m


#endif

