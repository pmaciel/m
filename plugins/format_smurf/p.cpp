
#include <cmath>
#include <cassert>
#include <iostream>
#include <vector>
#include <string>
#include "ext/GetPot"
#include "ext/fparser.hh"
#include "smurf.h"


namespace aux {

  // more helful zone structure, holding the variable values
  struct zone : public SmURF::TecZone {
    zone(const SmURF::TecZone& other) { this->operator=(other); }
    zone& operator=(const SmURF::TecZone& other) {this->SmURF::TecZone::operator=(other); return *this; }
    std::vector< std::vector< double > > vv;
  };

}


using namespace std;


int main(int argc, char **argv)
{
  // get options
  GetPot o(argc,argv);
  if (o.search(2,"-h","--help") || o.size()<=1) {
    cout << "binary tecplot ordered zones filter, usage: " << o[0] << " [options], with:" << endl
         << "  -h|--help             this help" << endl
         << endl
         << "  -i|--input-file       [inputfile] (default \"input.bin.plt\")" << endl
         << "  -o|--output-file      [outputfile] (default \"[inputfile].filter.bin.plt\")" << endl
         << "  -f|--filter           [function] filter to apply, default \"\")" << endl
         << endl
         << "  -of|--output-float    Tecplot data type (default double)" << endl
         << "  -or|--output-reverse  Tecplot reversed output (default direct)" << endl
         << "  -ov|--output-version  [number] Tecplot version (default 107)" << endl;
    return 0;
  }
  const string finp(o.follow("input.bin.plt",                 2,"-i","--input-file")),
               fout(o.follow((finp+".filter.bin.plt").c_str(),2,"-o","--output-file")),
               filt(o.follow("",                              2,"-f","--filter"));
  const SmURF::DataType outD = (o.search(    2,"-of","--output-float")? SmURF::FLOAT : SmURF::DOUBLE );
  const bool            outR =  o.search(    2,"-or","--output-reverse");
  const unsigned        outV =  o.follow(107,2,"-ov","--output-version");


  clog << "info: headers section: read..." << endl;
  string              tec_title;  // input file title
  vector< string >    tec_vars;   // input file variable names
  vector< aux::zone > tec_zones;  // input file zones descriptions
  SmURF::MeshReader sreader(finp);
  sreader.readMainHeader(tec_title,tec_vars);
  {
    vector< SmURF::TecZone > vtz = sreader.readZoneHeaders();
    tec_zones.assign(vtz.begin(),vtz.end());
  }
  clog << "info: headers section: read." << endl;


  clog << "info: setup filter..." << endl;
  FunctionParser filter;
  {
    filter.AddConstant("pi",M_PI);  // useful stuff!

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
  clog << "info: setup filter." << endl;


  clog << "info: data section: read and filter..." << endl;
  unsigned isum  = 0,  // total number of entries
           iskip = 0,  // skipped i entries
           zskip = 0;  // .. zone entries
  {
    vector< double > ev(tec_vars.size(),0.);  // i entries variable values
    vector< vector< unsigned > > _ve;         // matrix of elements nodes (dummy)
    for (vector< aux::zone >::iterator z=tec_zones.begin(); z!=tec_zones.end(); ++z) {
      sreader.readZoneData(*z,_ve,z->vv);
      isum  += z->i;
      iskip += z->i;
      if (z->type==SmURF::ORDERED && (z->i)>0) {

        // filter some entries out (in reverse, for performance)
        // (set entry values, evaluate condition and filter)
        for (int i=((int) z->i)-1; i>=0; --i) {
          for (unsigned v=0; v<tec_vars.size(); ++v)
            ev[v] = (z->vv)[v][i];
          if (!filter.Eval(&ev[0])) {
            for (unsigned v=0; v<tec_vars.size(); ++v)
              (z->vv)[v].erase((z->vv)[v].begin()+i);
            --z->i;
          }
        }

        iskip -= z->i;
      }
      zskip += (!z->i || z->type!=SmURF::ORDERED? 1:0);
    }
    cout << "info: filter"
         << " z:" << tec_zones.size() << '>' << tec_zones.size() - zskip
         << " i:" << isum             << '>' << isum-iskip << endl;
  }
  clog << "info: data section: read and filter." << endl;


  if (zskip>=tec_zones.size()) {
    clog << "warn: filter is never satisfied, skipping output file" << endl;
    return 1;
  }


  clog << "info: headers and data section: write..." << endl;
  SmURF::MeshWriter swriter(fout,outD,outR,outV);
  swriter.writeMainHeader(tec_title,tec_vars);

  for (vector< aux::zone >::iterator z=tec_zones.begin(); z!=tec_zones.end(); ++z)
    if (z->type==SmURF::ORDERED && z->i)
      swriter.writeZoneHeader(z->time,z->type,z->pack,z->title,z->i);

  vector< vector< unsigned > > _ve;  // matrix of elements nodes (dummy)
  for (vector< aux::zone >::iterator z=tec_zones.begin(); z!=tec_zones.end(); ++z)
    if (z->type==SmURF::ORDERED && z->i)
      swriter.writeZoneData(z->type,z->pack,_ve,z->vv,-1);
  clog << "info: headers and data section: write." << endl;


  return 0;
}
