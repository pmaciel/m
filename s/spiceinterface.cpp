
#include "mlog.h"
#include "spiceinterface.h"


namespace spice {


int interface_t::SendChar(char *msg, int id, void *data)
{
  const char *cerr, *cout, *m(
    (cerr=std::strstr(msg,"stderr "))!=NULL? cerr+7 :
    (cout=std::strstr(msg,"stdout "))!=NULL? cout+7 : msg );
  (cerr!=NULL? mlog::warn() : mlog::info()) << "spice[" << id << "]: " << m << mlog::nl;
  return 0;
}


int interface_t::SendStatus(char* msg, int id, void *data)
{
  mlog::info() << "spice[" << id << "]: " << msg << mlog::nl;
  return 0;
}


int interface_t::ControledExit(int st, bool unload, bool exit, int id, void *data)
{
  (st? mlog::warn() : mlog::info()) << "spice[" << id << "]:"
    << " status:" << st << " unload:" << unload << " exit:" << exit << mlog::nl;
  return 0;
}


int interface_t::BGThreadRunning(bool running, int id, void *data)
{
  mlog::info() << "spice[" << id << "]:" << " running:" << running << mlog::nl;
  return 0;
}


}  // namespace spice
