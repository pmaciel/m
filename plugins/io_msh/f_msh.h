#ifndef f_msh_h
#define f_msh_h

#include "mkernel.h"

namespace std { class fstream; }

// module for .msh
class f_msh : public m::mfinput {

 public:
  void read(GetPot& o, m::mmesh& m);

 private:
  void readMeshFormat          (const std::string& end, std::ifstream& f, m::mmesh& m);
  void readNodes               (const std::string& end, std::ifstream& f, m::mmesh& m);
  void readElements            (const std::string& end, std::ifstream& f, m::mmesh& m);
  void readPeriodic            (const std::string& end, std::ifstream& f, const m::mmesh& m);
  void readPhysicalNames       (const std::string& end, std::ifstream& f, m::mmesh& m);
  void readNodeData            (const std::string& end, std::ifstream& f, m::mmesh& m);
  void readElementData         (const std::string& end, std::ifstream& f, const m::mmesh& m);
  void readInterpolationScheme (const std::string& end, std::ifstream& f, const m::mmesh& m);
  void readComments            (const std::string& end, std::ifstream& f);
  void readError               (std::ifstream& f);

};

#endif

