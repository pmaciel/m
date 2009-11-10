
#include "mkernel.h"
#include "mfactory.h"


using namespace std;
using namespace m;


string extension(const string& fn)
{
  const string::size_type idx = fn.rfind('.');
  return fn.substr(idx!=string::npos? idx:0);
}


int main(int argc, char **argv)
{
  // options, formats and transformations
  GetPot o(argc,argv);
  mfactory< mfinput    >* fi = mfactory< mfinput    >::instance();
  mfactory< mfoutput   >* fo = mfactory< mfoutput   >::instance();
  mfactory< mtransform >* ft = mfactory< mtransform >::instance();


  // help (options and descriptions)
  if (o.search(2,"-h","--help") || o.size()<=1) {

    // set keys and descriptions
    vector< string > vk, vt;
    vk.push_back("");                    vt.push_back("");
    vk.push_back("General options:");    vt.push_back("");
    vk.push_back("-h|--help");           vt.push_back("this help");
    vk.push_back("-i [file]");           vt.push_back("input file(s)");
    vk.push_back("-o [file]");           vt.push_back("output file(s)");
    vk.push_back("-t[...]");             vt.push_back("transformations");
    vk.push_back("-version|--version");  vt.push_back("display version (-: proceed, --: stop)");
    vk.push_back("");                    vt.push_back("");
    vk.push_back("Input formats:");      vt.push_back("");
    fi->getkeys(vk);                     fi->gettxts(vt);
    vk.push_back("");                    vt.push_back("");
    vk.push_back("Output formats:");     vt.push_back("");
    fo->getkeys(vk);                     fo->gettxts(vt);
    vk.push_back("");                    vt.push_back("");
    vk.push_back("Transformations:");    vt.push_back("");
    ft->getkeys(vk);                     ft->gettxts(vt);

    // output with keys strings adjusted to the same length
    unsigned l = 0;
    for (unsigned i=0; i<vk.size(); ++i)
      l = (l>vk[i].length()? l:vk[i].length());
    l+=2;
    for (unsigned i=0; i<vk.size(); ++i)
      vk[i].insert(vk[i].end(),l-vk[i].length(),' ');

    cout << endl << "Mesh converter. Usage: " << o[0] << " [options]" << endl;
    for (unsigned i=0; i<vk.size(); ++i)
      cout << vk[i] << vt[i] << endl;
    exit(0);

  }
  else if (o.search(2,"-version","--version")) {

    cout << endl
         << "Mesh converter " << MVERSION << endl;

    cout << "mfinput keys:";
    vector< string > vk;
    fi->getkeys(vk);
    for (unsigned i=0; i<vk.size(); ++i)
      if (vk[i].length())
        cout << ' ' << vk[i];
    cout << endl;

    cout << "mfoutput keys:";
    vk.clear();
    fo->getkeys(vk);
    for (unsigned i=0; i<vk.size(); ++i)
      if (vk[i].length())
        cout << ' ' << vk[i];
    cout << endl;

    cout << "mtransform keys:";
    vk.clear();
    ft->getkeys(vk);
    for (unsigned i=0; i<vk.size(); ++i)
      if (vk[i].length())
        cout << ' ' << vk[i];
    cout << endl;

    if (o.search("--version"))
      exit(0);

  }


  // mesh structure
  mmesh m;


  // process command line arguments in order
  for (unsigned i=1; i<o.size(); ++i) {
    const string arg(o[i]);
    const string val(o[i+1]);
    o.set_cursor(i);
    if (arg=="-i") {

      // multiple input merges another mesh into current
      const string key(extension(val));
      if (fi->search(key.c_str())) {
        cout << "::read \"" << val << "\"..." << endl;
        mfinput* p = fi->Create(key);
        if (m.v()) {
          mmesh m2;
          p->read(o,m2);
          cout << "::merge [d/n/e]: " << m.d() << "+" << m2.d()
                             << " / " << m.n() << "+" << m2.n()
                             << " / " << m.e() << "+" << m2.e() << "..." << endl;
          m.merge(m2);
          cout << "::merge [d/n/e]: " << m.d()
                             << " / " << m.n()
                             << " / " << m.e() << "." << endl;
        }
        else {
          p->read(o,m);
        }
        delete p;
        cout << "::read \"" << val << "\"." << endl;
      }

    }
    else if (arg.find("-t")==0) {

      const string key(arg);
      if (ft->search(key.c_str())) {
        cout << "::transform \"" << key << "\"..." << endl;
        mtransform* p = ft->Create(key);
        p->transform(o,m);
        delete p;
        cout << "::transform \"" << key << "\"." << endl;
      }

    }
    else if (arg=="-o") {

      // write current mesh
      const string key(extension(val));
      if (fo->search(key.c_str())) {
        cout << "::write \"" << val << "\"..." << endl;
        mfoutput* p = fo->Create(key);
        p->write(o,m);
        delete p;
        cout << "::write \"" << val << "\"." << endl;
      }

    }
  }


  delete fi;
  delete fo;
  delete ft;


  return 0;
}

