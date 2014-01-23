
#include "mfactory.h"
#include "t_debug.h"

using namespace std;
using namespace m;


Register< mtransform,t_debug > mt_debug("-t42","intentionally undocumented feature :)");


void t_debug::transform(GetPot& o, mmesh& m, const XMLNode& x)
{
}

