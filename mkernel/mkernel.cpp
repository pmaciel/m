
#include <algorithm>  // for transform
#include <cctype>     // for toupper
#include <sstream>

#include "mkernel.h"


using namespace std;
using namespace m;


namespace m {
namespace utils {


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
  stringstream ss(s);
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


string trim(const string& s, const string& t)
{
  string str = s;
  return str.erase(s.find_last_not_of(t)+1)    // trim right
            .erase(0,s.find_first_not_of(t));  // trim left
}


string upper(const string& s)
{
  string r(s);
  std::transform(r.begin(),r.end(),r.begin(), (int(*)(int)) toupper);
  return r;
}


}  // namespace utils
}  // namespace m
