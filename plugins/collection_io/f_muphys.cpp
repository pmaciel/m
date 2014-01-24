
#include "mfactory.h"
#include "ext/xmlParser.h"
#include "f_muphys.h"

using namespace std;
using namespace m;


Register< mfoutput,f_muphys > mf_muphys(".muphys","MuPhyS output format");


void f_muphys::write(GetPot& o, const mmesh& m, const XMLNode& x)
{
  const string fn(string(o.get(o.inc_cursor(),""))+".xml");
  const unsigned d =  m.d();  // number of dimensions
  ostringstream s;            // utility stream

  // setup declaration and root node
  XMLNode xroot = XMLNode::createXMLTopNode("xml",TRUE),
          xmesh = xroot.addChild("spaceMesh");
  xroot.addAttribute("version","1.0");
  xroot.addAttribute("encoding","UTF-8");

  // coordinates
  s << d;
  xmesh.addChild("coordinates");
  xmesh.getChildNode("coordinates").addAttribute("type","Cartesian");
  xmesh.getChildNode("coordinates").addText(s.str().c_str());
  s.str("");

  // nodes
  xmesh.addChild("nodes");
  for (unsigned i=0; i<m.n(); ++i) {
    s << endl;
    for (unsigned j=0; j<d; ++j)
      s << ' ' << m.vv[j][i];
  }
  s << endl;
  xmesh.getChildNode("nodes").addText(s.str().c_str());
  s.str("");

  // zones
  xmesh.addChild("zoneMap");
  for (unsigned i=0; i<m.z(); ++i) {
    XMLNode z = xmesh.getChildNode("zoneMap").addChild("zone");
    z.addAttribute("label",m.vz[i].n.c_str());
    z.addAttribute("type",
      (m.vz[i].t==FELINESEG?     "Line_2"        :
      (m.vz[i].t==FETRIANGLE?    "Triangle_3"    :
      (m.vz[i].t==FETETRAHEDRON? "Tetrahedron_4" :
                                 "" ))));
    for (unsigned e=0; e<m.e(i); ++e) {
      const vector< unsigned >& en = m.vz[i].e2n[e].n;
      s << endl;
      for (vector< unsigned >::const_iterator n=en.begin(); n!=en.end(); ++n)
        s << ' ' << *n;
    }
    s << endl;
    z.addText(s.str().c_str());
    s.str("");
  }

  // write it
  xroot.writeToFile((fn).c_str());
  xroot.deleteNodeContent();
  xroot = XMLNode::emptyXMLNode;
}

