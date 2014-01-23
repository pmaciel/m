
#ifndef m_hal_h
#define m_hal_h


#include <string>
#include "ext/xmlParser.h"


// HAL: Hierarchical Abstraction Layer for generic data manipulation
namespace m {
  namespace hal {


// xobject: primitive type to wrap XMLNode and provide copy/assignment methods
class xobject : public XMLNode {

 public:

  // default contructor
  xobject(XMLCSTR tag="xobject") {
    XMLNode::operator=(XMLNode::createXMLTopNode(tag));
  }

  // interpret (recognized) attributes
  xobject& flatten();

  // append/update attributes (local only)
  xobject& operator+=(const XMLNode& x);

  // append/update attributes (nested attributes as well)
  xobject& operator*=(const XMLNode& x);

  // (XMLNode overrides)
  xobject(const XMLNode& x) { XMLNode::operator=(x); }
  xobject getChildNode(int i=0) const                   { return XMLNode::getChildNode(i);      }
  xobject getChildNode(XMLCSTR name, int i)  const      { return XMLNode::getChildNode(name,i); }
  xobject getChildNode(XMLCSTR name, int *i=NULL) const { return XMLNode::getChildNode(name,i); }

 private:

  /*
   * retrieve an xobject by path, following regular expression:
   * "^([^:]+):([^@]+)@([^=]+)=([^=]+)$", with
   * option 1: ($1)file:($2)tagp:($3)tagc@($4)alabel=($5)avalue
   * option 2: ($1)file:($2)tagp
   * option 3: ($1)file
   */
  static xobject getXObjectByPath(const std::string& path);

};


  }  // namespace hal
}  // namespace m


#endif

