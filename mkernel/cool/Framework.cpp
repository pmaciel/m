
#include <iostream>
#include "cool/Framework.hh"
using namespace std;


void cf_always_assert_desc(const string& m, const bool a)
{
  if (!a) {
    cerr << m << endl;
    throw 42;
  }
}


void cf_assert_desc(const string& m, const bool a)
{
  cf_always_assert_desc(m,a);
}


void cf_assert(const bool a)
{
  cf_assert_desc("cf_assert failed",a);
}

