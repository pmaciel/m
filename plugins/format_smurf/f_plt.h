#ifndef f_plt_h
#define f_plt_h

#include "mkernel.h"

struct TecZone {
  TecZone(
      const std::string _t="",
      const unsigned _nenodes=0,
      const unsigned _i=1,
      const unsigned _j=1,
      const m::mtype _et=m::ORDERED,
      const bool _isblock=false,
      const bool _isshared=false) :
    t(_t),
    nenodes(_nenodes),
    i(_i),
    j(_j),
    et(_et),
    isblock(_isblock),
    isshared(_isshared) {}
  std::string t;
  unsigned nenodes;
  unsigned i, j;
  m::mtype et;
  bool isblock;
  bool isshared;
};

// module for .plt
class f_plt : public m::mfinput,
              public m::mfoutput {
 public:
  void read(GetPot& o, m::mmesh& m);
  void write(GetPot& o, const m::mmesh& m);

 private:
  // PlatingMaster specific functions
  void in_platingmaster_2d(std::ifstream& f, m::mmesh& m);
  void in_platingmaster_3d(std::ifstream& f, m::mmesh& m);
  void getPMBoundaryZonePointNodeIndices(const m::mmesh& m, std::vector< unsigned >& vi, const std::vector< double >& vx, const std::vector< double >& vy);
  void getPMBoundaryZoneElementsFromNodeIndices(const m::mmesh& m, const std::vector< unsigned >& vi, const std::vector< double >& vx, const std::vector< double >& vy, std::vector< m::melem >& ve);

 public:
  // tecplot specific functions (input)
  static std::vector< std::string > getVariables(const std::string& s);
  static TecZone getZoneHeader(const std::string& s, std::string& n, m::mtype& t);
  static void readZoneNodeValues(std::ifstream& f, std::vector< std::vector< double > >& vv, unsigned N, unsigned Nvars, bool isblock);
  static void readZoneConnectivity(std::ifstream& f, std::vector< m::melem >& ve, unsigned N, unsigned Nenodes, m::mtype& t);

  // tecplot specific functions (output)
  static std::string setVariables(const std::vector< std::string >& vn);
  static std::string setZoneHeader(const TecZone& tz, const unsigned Nvars=0, const double solutiontime=0.);
  static void writeZoneNodeValues(std::ofstream& f, const std::vector< std::vector< double > >& vv, bool isblock, const int& sharezone=-1);
  static void writeZoneNodeValues(std::ofstream& f, const std::vector< double >& v);
  static void writeZoneConnectivity(std::ofstream& f, const std::vector< m::melem >& ve, const m::mtype& t);

 public:
  // string manipulation
  static std::string trimright(const std::string& s, const std::string& t=" ");
  static std::string trimleft(const std::string& s, const std::string& t=" ");
  static std::string trim(const std::string& s, const std::string& t=" ");
  static std::vector< std::string > splitstring(const std::string& s);
  static std::string upper(const std::string& s);
};

#endif

