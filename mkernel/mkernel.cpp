
#include <sstream>

#include "mkernel.h"


using namespace m;


namespace m {
namespace utils {


using std::pair;
using std::string;
using std::vector;


vector< pair< string,string > > getoperands(const string &s)
{
  //TODO make use of the (better) split methods

  // split string find a '=' then a ':'
  vector< pair< string,string > > r;
  string::size_type p1 = 0,
                         p2 = 0;
  while (p2!=string::npos) {
    pair< string,string > p("","");

    p2 = std::min< string::size_type >(s.find(":",p1),s.find("=",p1));
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


vector< string >& split(const string &s, char delim, vector< string > &elems)
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
