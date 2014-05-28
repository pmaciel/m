
#include <iostream>
#include <sstream>

#include "utils.h"


namespace m {
namespace utils {


using std::string;
using std::vector;


unsigned getvindex(const mmesh& m, const string& n)
{
  for (unsigned i=0; i<m.vn.size(); ++i)
    if (n==m.vn[i])
      return i;
  std::cerr << "error: variable name not found: \"" << n << "\"" << std::endl;
  throw 42;
  return 0;
}


unsigned getzindex(const mmesh& m, const string& n)
{
  for (unsigned i=0; i<m.vz.size(); ++i)
    if (n==m.vz[i].n)
      return i;
  std::cerr << "error: zone name not found: \"" << n << "\"" << std::endl;
  throw 42;
  return 0;
}


vector< std::pair< string, string > > getoperands(const string& s)
{
  //TODO make use of the (better) split methods, like in boost

  // split string find a '=' then a ':'
  vector< std::pair< string,string > > r;
  string::size_type p1 = 0,
                         p2 = 0;
  while (p2!=string::npos) {
    std::pair< string,string > p("","");

    p2 = std::min(s.find(":",p1),s.find("=",p1));
    p.first = s.substr(p1,(p2==string::npos? p2:p2-p1));
    if (s.find("=",p1)<s.find(":",p1)) {
      p1 = p2+1;
      p2 = s.find(":",p1);

      p.second = s.substr(p1,(p2==string::npos? p2:p2-p1));
      /* older version, maybe works better
        string s2 = s.substr(p1,(p2==string::npos? p2:p2-p1));
        istringstream ss(s2);
        ss >> p.second;
      */
    }

    p1 = p2+1;
    r.push_back(p);
  }
  return r;
}


vector< string >& split(const string &s, char delim, vector< string >& elems)
{
  std::stringstream ss(s);
  string item;
  while (std::getline(ss,item,delim))
    elems.push_back(item);
  return elems;
}


vector< string > split(const string &s, char delim)
{
  vector< string > elems;
  split(s,delim,elems);
  return elems;
}


}  // namespace utils
}  // namespace m
