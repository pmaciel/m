
#include <cstdio>  // for sscanf
#include <iostream>
#include "hal.h"


namespace m {
  namespace hal {


#if 0
// hsh instance
hsh* hsh::m_instance;
#endif


xobject& xobject::flatten()
{
  // interpret xobject attribute:
  // 1. get & delete attribute, to avoid conflicts when nesting
  // 2. get new node, swap with internal (replace) and add nested properties
  while (isAttributeSet("xobject")) {

    const std::string str(getAttribute("xobject"));
    deleteAttribute("xobject");
    if (!str.length())
      continue;

    xobject y(getXObjectByPath(str));
    if (!y.isEmpty()) {
      std::swap(*this,y);
      operator*=(y);
    }

  }

  // interpret nested attributes (recursive)
  for (int i=0; i<nChildNode(); ++i)
    xobject::getChildNode(i).flatten();

  // return itself
  return *this;
}


xobject& xobject::operator+=(const XMLNode& x)
{
  for (int i=0; i<x.nAttribute(); ++i)
    updateAttribute(x.getAttributeValue(i),NULL,x.getAttributeName(i));
  return *this;
}


xobject& xobject::operator*=(const XMLNode& x)
{
  //TODO: allow matching with empty nodes, would be a much better way
  // (by path, perhaps?)

#if 1
  // xml node structure must match before continuing
  bool match(
    std::string(getName())==x.getName() &&
    nChildNode()==x.nChildNode() );
  for (int i=0; i<x.nChildNode() && match; ++i)
    match = std::string(getChildNode(i).getName())==x.getChildNode(i).getName();
  if (!match) {
    std::cerr << "xobject: xml node structure doesn't match!" << std::endl;
    throw 42;
  }

  // update the attributes, then recurse into child nodes
  operator+=(x);
  for (int i=0; i<x.nChildNode(); ++i)
    getChildNode(i) *= x.getChildNode(i);
#else
  // update the attributes, then set child nodes attributes
  operator+=(x);
  for (int i=0; i<x.nElement(); ++i) {
    XMLNodeContents c = x.enumContents(i);
    if (c.etype==eNodeChild) {
      XMLNode xpto = c.child;
      XMLSTR str = xpto.createXMLString(0);
      std::cout << "child:\"" << str << '"' << std::endl;
      freeXMLString(str);
    }
  }
#endif

  // return itself
  return *this;
}


xobject xobject::getXObjectByPath(const std::string& path)
{
  char r[5][50];
  int rc;
  if      ((rc=sscanf(path.c_str(),"%50[^:]:%50[^:]:%50[^@]@%50[^=]=%50[^=]",r[0],r[1],r[2],r[3],r[4]))==5) {}
  else if ((rc=sscanf(path.c_str(),"%50[^:]:%50[^:]",                        r[0],r[1]               ))==2) {}
  else if ((rc=sscanf(path.c_str(),"%50[^:]",                                r[0]                    ))==1) {}
  else { rc=0; }
  const std::string file  (rc>0? r[0]:""),
                    tagp  (rc>1? r[1]:""),
                    tagc  (rc>2? r[2]:""),
                    alabel(rc>3? r[3]:""),
                    avalue(rc>4? r[4]:"");
  xobject y(file.length()? XMLNode::openFileHelper(file.c_str())
                         : XMLNode::emptyNode() );
  y = (y.isEmpty()?    y :
      (tagc.length() && alabel.length() && avalue.length()? y.getChildNodeByPath(tagp.c_str()).getChildNodeWithAttribute(tagc.c_str(),alabel.c_str(),avalue.c_str()) :
      (tagp.length()?  y.getChildNodeByPath(tagp.c_str()) :
      (y.nChildNode()? y.getChildNode(0) :
                       XMLNode::emptyNode() ))));
  if (y.isEmpty()) {
    std::cerr << "xobject: empty xml node, likely problems with xobject=\"\" syntax!" << std::endl;
    throw 42;
  }
  return y;
}


  }  // namespace hal
}  // namespace m

