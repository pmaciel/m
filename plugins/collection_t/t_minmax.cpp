
#include <sstream>
#include <algorithm>
#include "mfactory.h"
#include "t_minmax.h"

using namespace std;
using namespace m;


Register< mtransform,t_minmax > mt_minmax("-tminmax","[str:...] [str:...] check variable minimum/maximum values at given zones");


void t_minmax::transform(GetPot& o, mmesh& m, const XMLNode& x)
{
  cout << "::minmax..." << endl;

  // find variables/zones to check
  // (a parameter of ":" makes it check all variables/zones)
  string vstr = o.get(o.inc_cursor(),""),
         zstr = o.get(o.inc_cursor(),"");
  vector< bool > chkv(m.v(),vstr==":"),
                 chkz(m.z(),zstr==":");
  if (vstr!=":") {
    replace(vstr.begin(),vstr.end(),':',' ');
    istringstream is(vstr);
    string name;
    while (is >> name)
      for (unsigned j=0; j<m.v(); ++j)
        chkv[j] = chkv[j] || m.vn[j]==name;
  }
  if (zstr!=":") {
    replace(zstr.begin(),zstr.end(),':',' ');
    istringstream is(zstr);
    string name;
    while (is >> name)
      for (unsigned j=0; j<m.z(); ++j)
        chkz[j] = chkz[j] || m.vz[j].n==name;
  }


  if (!m.z()) {

    // if no zones are present, check all entries in point cloud
    for (unsigned j=0; j<m.v(); ++j) {
      if (!chkv[j])
        continue;
      const vector< double >& v = m.vv[j];  // variable shortcut
      cout << "  variable \"" << m.vn[j] << "\""
           << "  min/max: " << *min_element(v.begin(),v.end())
                   << " / " << *max_element(v.begin(),v.end()) << endl;
    }

  }
  else {

    // check by zone, then by variable
    vector< bool > chknode(m.n());  // nodes to check
    for (unsigned i=0; i<m.z(); ++i) {
      if (!chkz[i]) continue;
      const mzone& z = m.vz[i];  // zone shortcut

      // mark nodes to check
      chknode.assign(m.n(),false);
      for (vector< melem >::const_iterator e=z.e2n.begin(); e!=z.e2n.end(); ++e)
        for (vector< unsigned >::const_iterator n=(e->n).begin(); n!=(e->n).end(); ++n)
          chknode[*n] = true;
      if (!count(chknode.begin(),chknode.end(),true)) {
        cout << "::minmax: warning: zone \""     << z.n     << "\" doesn't have any nodes, skipping" << endl;
        continue;
      }

      for (unsigned j=0; j<m.v(); ++j) {
        if (!chkv[j]) continue;

        // check for minimum/maximum
        double vmin =  1.7e308,
               vmax = -1.7e308;
        const vector< double >& v = m.vv[j];  // variable shortcut
        for (unsigned k=0; k<m.n(); ++k)
          if (chknode[k]) {
            vmin = v[k]<vmin? v[k]:vmin;
            vmax = v[k]>vmax? v[k]:vmax;
          }

        // summary
        cout << "  zone \""     << z.n     << "\""
             << "  variable \"" << m.vn[j] << "\""
             << "  min/max: " << vmin << " / " << vmax << endl;
      }
    }
  }

  cout << "::minmax." << endl;
}

