
#include <sstream>
#include "mlog.h"
#include "spiceinterface.h"


namespace spice {


int interface_t::SendChar(char *msg, int id, void *data)
{
  std::istringstream is(msg);
  std::string tag;
  is >> tag;
  (tag=="stderr"? mlog::warn() : mlog::info()) << "spice[" << id << "]: "
    << (tag=="stderr" || tag=="stdout"? is.str().substr(tag.length()) : msg)
    << mlog::nl;
  return 0;
}


int interface_t::SendStatus(char* msg, int id, void *data)
{
  mlog::info() << "spice[" << id << "]: " << msg << mlog::nl;
  return 0;
}


int interface_t::ControledExit(int st, bool dounload, bool doexit, int id, void *data)
{
  std::ostringstream is;
  st? is << " status:" << st : is;
  (st? mlog::warn() : mlog::info()) << "spice[" << id << "]:" << is.str()
    << " unload:" << (dounload? "true":"false")
    << " exit:"   << (doexit?   "true":"false")
    << mlog::nl;
  return 0;
}


int interface_t::BGThreadRunning(bool isrunning, int id, void *data) {
  mlog::info() << "spice[" << id << "]:"
    << " running:" << (isrunning? "true":"false")
    << mlog::nl;
  return 0;
}


}  // namespace spice
