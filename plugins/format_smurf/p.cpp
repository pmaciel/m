
#include <cmath>
#include <cassert>
#include <iostream>
#include <vector>
#include <string>
#include "ext/GetPot"
#include "ext/fparser.hh"
#include "smurf.h"


using namespace std;


int main(int argc, char **argv)
{
  // get options
  GetPot o(argc,argv);
  if (o.search(2,"-h","--help") || o.size()<=1) {
    cout << "smurf filter, usage: " << o[0] << " [options], with:" << endl
         << "  -h|--help                       this help" << endl
         << "  -if|--input-file [inputfile]    (default \"input.smurf\")" << endl
         << "  -of|--output-file [outputfile]  (default \"[inputfile].filter.smurf\")" << endl
         << "  -f|--filter [str]               filter to apply (function of file variables, default \"\")" << endl;
    return 0;
  }
  const string finp(o.follow("input.smurf",2,"-if","--input-file")),
               fout(o.follow((finp+".filter.smurf").c_str(),2,"-of","--output-file")),
               filt(o.follow("",2,"-f","--filter"));
  if (!filt.size()) {
    cerr << "error: no filter, nothing to do!" << endl;
    return 1;
  }


  // setup reader and writer
  SmURF::MeshReader sreader(finp);
  SmURF::MeshWriter swriter(fout);


  // setup function, by setting variable names and filter
  FunctionParser filter;
  filter.AddConstant("pi",M_PI);  // useful stuff!
  string           tec_title;     // input file title
  vector< string > tec_vars;      // input file variable names


  // main header section (read and write)
  sreader.readMainHeader(tec_title,tec_vars);
  swriter.writeMainHeader(tec_title,tec_vars);
  {
    assert(tec_vars.size());
    string vars(tec_vars[0]);  // variable names as a single string
    for (int i=1; i<(int) tec_vars.size(); ++i)
      vars += ',' + tec_vars[i];

    clog << "info: variables: \"" << vars << '"' << endl;
    clog << "info: function: \"" << filt << '"' << endl;
    const int r = filter.Parse(filt,vars);
    if (r>=0) {
      cerr << string(r+17,' ') << "^ " << filter.ErrorMsg() << endl;
      return 1;
    }
  }


  // zone headers section (read)
  vector< SmURF::TecZone > inp_zheaders = sreader.readZoneHeaders();


  // set data to store to write later
  vector< SmURF::TecZone >             out_zheaders;
  vector< vector< vector< double > > > out_vvv;


  // zone data section (read and filter)
  {
    vector< double > entry_vars(tec_vars.size(),0.);
    vector< vector< unsigned > > _ve;  // matrices of elements nodes
    vector< vector< double   > > _vv;  // matrices of values
    for (vector< SmURF::TecZone >::iterator z=inp_zheaders.begin(); z!=inp_zheaders.end(); ++z) {
      sreader.readZoneData(*z,_ve,_vv);
      if (z->type==SmURF::ORDERED) {

        // give this zone some room
        out_zheaders.push_back(*z);
        out_vvv.push_back(_vv);
        SmURF::TecZone&             filt_z  = out_zheaders.back();
        vector< vector< double > >& filt_vv = out_vvv.back();

        // filter some entries out (in reverse, for performance)
        for (unsigned i=filt_z.i; i>=0; --i) {
          for (unsigned v=0; v<tec_vars.size(); ++v)
            entry_vars[v] = filt_vv[v][i];
          if (!filter.Eval(&entry_vars[0])) {
            for (unsigned v=0; v<tec_vars.size(); ++v)
              filt_vv[v].erase(filt_vv[v].begin()+i);
            --filt_z.i;
          }
        }

      }
    }
  }


  // zone headers section (write)
  for (vector< SmURF::TecZone >::iterator z=out_zheaders.begin(); z!=out_zheaders.end(); ++z)
    swriter.writeZoneHeader(SmURF::ORDERED,SmURF::BLOCK,z->title,z->i);


  // zone data section (write)
  {
    vector< vector< unsigned > > _ve;  // matrices of elements nodes
    for (unsigned i=0; i<out_zheaders.size(); ++i)
      swriter.writeZoneData(out_zheaders[i].type,SmURF::BLOCK,_ve,out_vvv[i],-1);
  }


  return 0;
}


#if 0
void f_smurf::write(GetPot& o, const mmesh& m)
{
  const string fn(o.get(o.inc_cursor(),""));

}
#endif

