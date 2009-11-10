
#ifndef mfactory_h
#define mfactory_h

#include <iostream>
#include <vector>
#include <cstdarg>


namespace m {


/**
 * this is a factory (with a singleton pattern) made simple. you create one
 * "factory" (of class X) -- by accessing its instance -- and register multiple
 * "production lines" (of subclasses of X), then create "products" given a key.
 *
 * destruction of "production lines" is also properly handled.
 */

template< class TheProductType >
class AProductionLine {
 public:
  AProductionLine(const std::string& _key, const std::string& _txt) : key(_key), txt(_txt) {}
  virtual ~AProductionLine() {}
  virtual TheProductType* create() = 0;
  std::string key;
  std::string txt;
};


template< class TheProductType, class TheProduct >
class ProductionLine : public AProductionLine< TheProductType > {
 public:
  ProductionLine(const std::string& _key, const std::string& _txt) : AProductionLine< TheProductType >(_key,_txt) {}
  TheProductType* create() { return new TheProduct(); }
};


template< class TheProductType >
class mfactory {

 private:

  // singleton instance
  static mfactory< TheProductType >* m_instance;


 protected:

  // constructor for singleton
  mfactory() {}


 public:

  // access instance
  static mfactory< TheProductType >* instance() {
    if (m_instance==0)
      m_instance = new mfactory< TheProductType >();
    return m_instance;
  }


  // get keys
  void getkeys(std::vector< std::string >& v) {
    for (unsigned i=0; i<m_vlines.size(); ++i)
      v.push_back(m_vlines[i]->key);
  }


  // get keys descriptions
  void gettxts(std::vector< std::string >& v) {
    for (unsigned i=0; i<m_vlines.size(); ++i)
      v.push_back(m_vlines[i]->txt);
  }


  // destroy all production lines
  ~mfactory() {
    for (unsigned i=0; i<m_vlines.size(); ++i)
      delete m_vlines[i];
  }


  // search for a key
  bool search(const char* k) {
    const std::string key(k);
    if (!key.length())
      return false;
    for (unsigned i=0; i<m_vlines.size(); ++i)
      if (m_vlines[i]->key==k)
        return true;
    return false;
  }


  // search for one of multiple keys
  bool search(const unsigned n, const char* k, ...) {
    // search on empty list and first argument
    if (!n)         return false;
    if (search(k))  return true;

    // interpret variable argument list
    bool yes = false;
    va_list ap;
    va_start(ap,k);
    for(unsigned i=1; i<n && !yes; ++i) {
      const char* k = va_arg(ap,char*);
      yes = search(k);
    }
    va_end(ap);
    return yes;
  }


  // register key
  bool Register(AProductionLine< TheProductType >* l) {
    const std::string key = l->key;
    if (search(key.c_str())) {
      std::cout << "Warning: \"" << key << "\" already exists!" << std::endl;
      return false;
    }
    m_vlines.push_back(l);
    return true;
  }


  // unregister key
  void Unregister(const std::string& key) {
    for (unsigned i=0; i<m_vlines.size(); ++i) {
      if (m_vlines[i]->key==key) {
        delete m_vlines[i];
        m_vlines.erase(i);
        return;
      }
    }
    std::cout << "Warning: \"" << key << "\" unregistration not possible!" << std::endl;
  }


  // create product from key
  TheProductType* Create(const std::string& key) {
    for (unsigned i=0; i<m_vlines.size(); ++i)
      if (key.length() && m_vlines[i]->key==key)
        return m_vlines[i]->create();
    std::cout << "Warning: \"" << key << "\" creation not possible!" << std::endl;
    return NULL;
  }


 private:

  // production lines
  std::vector< AProductionLine< TheProductType >* > m_vlines;

};


/**
 * factory singleton instance
 */
template< class TheProductType >
mfactory< TheProductType >* mfactory< TheProductType >::m_instance;


/**
 * factory registration by object construction (to simplify factory use)
 */
template< class TheProductType, class TheProduct >
class Register {
 public:

  // register one key
  Register(const char* k, const char* t) {
    Register(1,k,t);
  }

  // register multiple keys
  Register(const unsigned n, const char* k, const char* t, ...) {

    // get variable arguments as (k,t) pairs
    std::vector< std::string > vk,
                               vt;
    vk.push_back(k);
    vt.push_back(t);
    va_list ap;
    va_start(ap,t);
    for(unsigned i=1; i<n; ++i) {
      const char* k = va_arg(ap,char*);  vk.push_back(k);
      const char* t = va_arg(ap,char*);  vt.push_back(t);
    }
    va_end(ap);

    // register (k,t) pairs
    mfactory< TheProductType >* f = mfactory< TheProductType >::instance();
    for(unsigned i=0; i<vk.size(); ++i)
      f->Register(new ProductionLine< TheProductType,TheProduct >(vk[i],vt[i]));
  }

};


/**
 * factory production by static function call (to simplify factory use)
 */
template< class TheProductType >
static TheProductType* Create(const std::string& key) {
  mfactory< TheProductType >* f = mfactory< TheProductType >::instance();
  return f->Create(key);
}


}  // namespace m

#endif

