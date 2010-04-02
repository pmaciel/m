
#include "mfactory.h"
#include "ext/xmlParser.h"
#include "f_muphys.h"

using namespace std;
using namespace m;


Register< mfoutput,f_muphys > mf_muphys(".muphys","MuPhyS output format");


void f_muphys::write(GetPot& o, const mmesh& m)
{
  const string fn(o.get(o.inc_cursor(),""));
  const unsigned d =  m.d();  // number of dimensions
  ostringstream s;            // utility stream

  // setup declaration and root node
  XMLNode r = XMLNode::createXMLTopNode("xml",TRUE),
          x = r.addChild("spaceMesh");
  r.addAttribute("version","1.0");
  r.addAttribute("encoding","UTF-8");

  // coordinates
  s << d;
  x.addChild("coordinates");
  x.getChildNode("coordinates").addAttribute("type","Cartesian");
  x.getChildNode("coordinates").addText(s.str().c_str());
  s.str("");

  // nodes
  x.addChild("nodes");
  for (unsigned i=0; i<m.n(); ++i) {
    s << endl;
    for (unsigned j=0; j<d; ++j)
      s << ' ' << m.vv[j][i];
  }
  s << endl;
  x.getChildNode("nodes").addText(s.str().c_str());
  s.str("");

  // zones
  x.addChild("zoneMap");
  for (unsigned i=0; i<m.z(); ++i) {
    XMLNode z = x.getChildNode("zoneMap").addChild("zone");
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
  r.writeToFile((fn+".xml").c_str());
  r.deleteNodeContent();
  r = XMLNode::emptyXMLNode;
}

