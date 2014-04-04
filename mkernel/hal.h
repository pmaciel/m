
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


#if 0
// xpath: primitive type for path properties, describing an hierarchy entry
struct xpath : xobject {
  xpath(XMLCSTR s) : xobject() {
    updateName("hal::xpath");
    updateAttribute(s,NULL,"hal::name");
    updateAttribute("",NULL,"hal::type");
  }
  std::string name(const std::string& s="") { if (s.length()) updateAttribute(s.c_str(),NULL,"hal::name"); return getAttribute< std::string >("hal::name"); }
  std::string type(const std::string& s="") { if (s.length()) updateAttribute(s.c_str(),NULL,"hal::type"); return getAttribute< std::string >("hal::type"); }
/*
  template< typename T > const T&  r() { return w< T >(); }
  template< typename T >       T&  w() { return *(new T); }
                               int x() { return 0; }
*/
};


// HAL node type describing arbitrary data and its properties
struct node : xpath
{
  node(XMLCSTR s) : xpath(s) { updateName("hal::node"); }
};


// HAL link type describing a link to directories and data
struct link : xpath
{
  link(XMLCSTR s) : xpath(s) { updateName("hal::link"); }
};


// HAL directory/category type to hold other HAL types
struct directory : xpath
{
  directory(XMLCSTR s) : xpath(s) { updateName("hal::directory"); }
};


// HAL shell providing an environment context and access to HAL entries
class hsh {

 public:

  // access instance
  static hsh* instance() { return (m_instance? m_instance : new hsh()); }

  // create node/directory
  int touch(XMLCSTR s) {
    if (!m_root.getChildNode(s).isEmpty())
      m_root.addChild(node(s));
    return 0;
  }
  int mkdir(XMLCSTR s) {
    if (m_root.getChildNode(s).isEmpty())
      return -1;
    m_root.addChild(directory(s));
    return 0;
  }
  int cd(XMLCSTR s="~") {
    if      (std::string(s)=="~")  {}
    else if (std::string(s)==".")  {}
    else if (std::string(s)=="..") {}
    return 0;
  }

  XMLCSTR getenv(XMLCSTR variable) {
    return m_environment.getAttribute(variable);
  }

 private:

  // singleton instance and constructor
  static hsh* m_instance;
  hsh() : m_root("/") {
    m_environment.updateAttribute("/",NULL,"PWD");
    mkdir("home");
    mkdir("lib");
    mkdir("bin");
    cd("home");
  }

public:
  // environment and root node
  xobject   m_environment;
  directory m_root;

};
#endif


  }  // namespace hal
}  // namespace m


#endif

