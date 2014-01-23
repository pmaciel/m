#ifndef t_meshdef_h
#define t_meshdef_h

#include "mkernel.h"


// module with a 2D Laplace equation solver
class t_meshdef : public m::mtransform {
 public:
  void transform(GetPot& o, m::mmesh& m, const XMLNode& x);
 private:
  // mesh utilities
//  std::vector< double > getvvalues(const std::string& s);

  template< typename T >
  std::vector< T > getvvalues(const std::string& s)
  {
    std::vector< std::string > q;
    boost::algorithm::split(q,s,boost::is_any_of(":"),boost::token_compress_off );
    std::vector< T > r;
    BOOST_FOREACH(const std::string& str,q) {
      r.push_back(boost::lexical_cast< T >(str.c_str()));
    }
    return r;
  }

  std::vector< m::mzone >::const_iterator getzoneit  (const std::string& n, const m::mmesh& m);
  unsigned getzoneidx (const std::string& n, const m::mmesh& m);
  unsigned getvaridx  (const std::string& n, const m::mmesh& m);
};

#endif

