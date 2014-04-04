
#include <iostream>
#include <vector>
#include <algorithm>
#include "ls_samg_parm.h"


namespace m {
  namespace aux {


void s_parm::set(const std::string &_s)
{
  // parameters set by functions
  if      (fpi!=NULL) { int    v = cvt< int    >(_s); (*fpi)( &v ); }
  else if (fpd!=NULL) { double v = cvt< double >(_s); (*fpd)( &v ); }
  else if (fps!=NULL) {
    int w = (int) _s.length();  // convert string to vector of int's...
    std::vector< int > v(w);
    for (int j=0; j<w; ++j)
      v[j] = (int) _s[j];
    (*fps)( &v[0],&w );        // ... then set the option
  }

  // parameters set "directly"
  else if (pd!=NULL)        { *pd = cvt< double >(_s); }
  else if (pi!=NULL && !p1) { *pi = cvt< int    >(_s); }

  // integer sub-parameters
  else if (pi!=NULL) {

    // build vectors of the sub-parameter and parameter numbers digits
    const int sign = (*pi<0? -1:1);
    std::vector< unsigned > parm,
                            subparm;
    {
      unsigned v = cvt< unsigned >(_s);
      while (v) {
        subparm.push_back(v % radix);
        v /= radix;
      }
      while (p2-p1+1>(int) subparm.size())
        subparm.push_back(0);
      std::reverse(subparm.begin(),subparm.end());
      if (p2 && (int) subparm.size()!=p2-p1+1) {
        std::cerr << "warning: sub-parameter \"" << n << "\" should have " << (p2-p1+1) << " digit(s)!" << std::endl;
        return;
      }

      unsigned numberp = std::max(*pi,-(*pi));
      while (numberp) {
        parm.push_back(numberp % radix);
        numberp /= radix;
      }
      std::reverse(parm.begin(),parm.end());
      while ((int) parm.size()<(int) subparm.size()+(p1-1))
        parm.push_back(0);
    }

    // set sub-parameter in the vector of digits
    for (int j=0; j<(int) subparm.size(); ++j)
      parm[p1+j-1] = subparm[j];

    // rebuild number from vector of digits
    int number = 0,
        base   = sign;
    while (parm.size()) {
      number += parm.back()*base;
      parm.pop_back();
      base *= radix;
    }
    *pi = number;

  }
}


void s_parm::def()
{
  if      (pd!=NULL)         *pd = d;
  else if (pi!=NULL && !p1)  *pi = i;
}


std::string s_parm::get()
{
  if      (pd!=NULL) { std::ostringstream oss;  oss << *pd;  return oss.str(); }
  else if (pi!=NULL) { std::ostringstream oss;  oss << *pi;  return oss.str(); }
  else return std::string();
}


  }  // namespace aux
}  // namespace m

