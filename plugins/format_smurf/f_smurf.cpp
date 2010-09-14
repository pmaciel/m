
#include <sstream>
#include "smurf.h"
#include "mfactory.h"
#include "f_smurf.h"

using namespace std;
using namespace m;


Register< mfinput,f_smurf > mf_smurf1(".smurf","Tecplot binary input format");
Register< mfoutput,f_smurf > mf_smurf2(6,".smurf","Tecplot binary output format",
                                         "","--smurf-point: output in point format (default: block)",
                                         "","--smurf-float: output in float precision (default: double)",
                                         "","--smurf-reverse: output reversed bytes (default: no)",
                                         "","--smurf-version [int]: output version (default: 111)",
                                         "","--smurf-solutiontime [real]: solution time (default: 0.)");


void f_smurf::read(GetPot& o, mmesh& m)
{
  // reader
  const string fn(o.get(o.inc_cursor(),""));
  SmURF::MeshReader mreader(fn);

  // headers section
  string title;
  mreader.readMainHeader(title,m.vn);
  vector< SmURF::TecZone > zheaders = mreader.readZoneHeaders();

  // data section
  m.vz.clear();
  for (vector< SmURF::TecZone >::iterator zone=zheaders.begin(); zone!=zheaders.end(); ++zone) {
    cout << "info: ZONE TITLE=\"" << zone->title << "\""
         << " ZONETYPE=" << (zone->type==SmURF::ORDERED?         "ORDERED":
                            (zone->type==SmURF::FELINESEG?       "FELINESEG":
                            (zone->type==SmURF::FETRIANGLE?      "FETRIANGLE":
                            (zone->type==SmURF::FEQUADRILATERAL? "FEQUADRILATERAL":
                            (zone->type==SmURF::FETETRAHEDRON?   "FETETRAHEDRON":
                            (zone->type==SmURF::FEBRICK?         "FEBRICK":"error") )))))
         << " DATAPACKING=" << (zone->pack==SmURF::BLOCK? "BLOCK":
                               (zone->pack==SmURF::POINT? "POINT":"error") )
         << " I,J,K=" << zone->i << "," << zone->j << "," << zone->k << endl;

    mzone z;
    z.n = zone->title;
    switch (zone->type) {
      case SmURF::ORDERED:          { z.t=ORDERED;         } break;
      case SmURF::FELINESEG:        { z.t=FELINESEG;       } break;
      case SmURF::FETRIANGLE:       { z.t=FETRIANGLE;      } break;
      case SmURF::FEQUADRILATERAL:  { z.t=FEQUADRILATERAL; } break;
      case SmURF::FETETRAHEDRON:    { z.t=FETETRAHEDRON;   } break;
      case SmURF::FEBRICK:          { z.t=FEBRICK;         } break;
      default:                      { z.t=ORDERED;         } break;
    }
    m.vz.push_back(z);

    vector< vector< unsigned > > ve;
    vector< vector< double   > > vv;
    mreader.readZoneData(*zone,ve,vv);

    // correct connectivity numbering
    if (ve.size() && z.t==FEBRICK) {
      if (ve[0][2]==ve[0][3] && ve[0][6]==ve[0][7]) { //*2
        m.vz.back().t = PRISM3;
        for (unsigned c=0; c<ve.size(); ++c) {
          const vector< unsigned >  ent = ve[c];  // element nodes, original (tecplot)
                vector< unsigned >& eng = ve[c];  // ...,           modified (gambit)
          eng.resize(6);
          eng[0] = ent[0];
          eng[1] = ent[1];
          eng[2] = ent[2];
          eng[3] = ent[4];
          eng[4] = ent[5];
          eng[5] = ent[6];
        }
      }
      else if ((ve[0][4]==ve[0][5]) && (ve[0][4]==ve[0][6]) && (ve[0][4]==ve[0][7])) { //*3
        m.vz.back().t = PYRAMID4;
        for (unsigned c=0; c<ve.size(); ++c) {
          const vector< unsigned >  ent = ve[c];  // element nodes, original (tecplot)
                vector< unsigned >& eng = ve[c];  // ...,           modified (gambit)
          eng.resize(5);
          eng[0] = ent[0];
          eng[1] = ent[1];
          eng[2] = ent[3];
          eng[3] = ent[2];
          eng[4] = ent[4];
        }
      }
      else { //*1
        for (unsigned c=0; c<ve.size(); ++c) {
          swap(ve[c][2],ve[c][3]);
          swap(ve[c][6],ve[c][7]);
        }
      }
    }
    m.vz.back().e2n = convert_to_vtelem(ve);

    if (m.vz.size()==1)  // first zone has the variable data
      m.vv.swap(vv);
  }
}


void f_smurf::write(GetPot& o, const mmesh& m)
{
  const string fn(o.get(o.inc_cursor(),""));

  // options & writer
  const SmURF::ZonePack pack     = (o.search("--smurf-point")? SmURF::POINT : SmURF::BLOCK );
  const SmURF::DataType datatype = (o.search("--smurf-float")? SmURF::FLOAT : SmURF::DOUBLE );
  const bool     reverse      = o.search("--smurf-reverse");
  const unsigned version      = o.follow(111, "--smurf-version");
  const double   solutiontime = o.follow(0.,  "--smurf-solutiontime");
  SmURF::MeshWriter mwriter(fn,datatype,reverse,version);

  // headers section
  mwriter.writeMainHeader("untitled",m.vn);
  for (unsigned i=0; i<m.z(); ++i) {
    const mzone& z = m.vz[i];
    const SmURF::ZoneType type = (z.t==ORDERED?         SmURF::ORDERED         :
                                 (z.t==FELINESEG?       SmURF::FELINESEG       :
                                 (z.t==FETRIANGLE?      SmURF::FETRIANGLE      :
                                 (z.t==FEQUADRILATERAL? SmURF::FEQUADRILATERAL :
                                 (z.t==FETETRAHEDRON?   SmURF::FETETRAHEDRON   :
                                 (z.t==FEBRICK?         SmURF::FEBRICK         : //*1
                                 (z.t==FEPOLYGON?       SmURF::FEPOLYGON       :
                                 (z.t==FEPOLYHEDRON?    SmURF::FEPOLYHEDRON    :
                                 (z.t==PRISM3?          SmURF::FEBRICK         : //*2
                                 (z.t==PYRAMID4?        SmURF::FEBRICK         : //*3
                                                        SmURF::ORDERED ))))))))));
    mwriter.writeZoneHeader(solutiontime,type,pack,z.n,m.n(),m.e(i));
  }
  if (!m.z() && m.v()) {
    // no connectivities present, but there is a point cloud
    mwriter.writeZoneHeader(solutiontime,SmURF::ORDERED,pack,"point_cloud",m.n());
  }

  // data section
  for (unsigned i=0; i<m.z(); ++i) {
    const mzone& z = m.vz[i];
    const SmURF::ZoneType type = (z.t==ORDERED?         SmURF::ORDERED         :
                                 (z.t==FELINESEG?       SmURF::FELINESEG       :
                                 (z.t==FETRIANGLE?      SmURF::FETRIANGLE      :
                                 (z.t==FEQUADRILATERAL? SmURF::FEQUADRILATERAL :
                                 (z.t==FETETRAHEDRON?   SmURF::FETETRAHEDRON   :
                                 (z.t==FEBRICK?         SmURF::FEBRICK         : //*1
                                 (z.t==FEPOLYGON?       SmURF::FEPOLYGON       :
                                 (z.t==FEPOLYHEDRON?    SmURF::FEPOLYHEDRON    :
                                 (z.t==PRISM3?          SmURF::FEBRICK         : //*2
                                 (z.t==PYRAMID4?        SmURF::FEBRICK         : //*3
                                                        SmURF::ORDERED ))))))))));

    // correct connectivity numbering
    vector< vector< unsigned > > ve = convert_from_vtelem(z.e2n);
    if (z.t==PRISM3) { //*2
      // these are represented by a FEBRICK, with coalesced nodes
      for (unsigned c=0; c<ve.size(); ++c) {
        const vector< unsigned >  eng = ve[c];  // element nodes, original (gambit)
              vector< unsigned >& ent = ve[c];  // ...,           modified (tecplot)
        ent.resize(8);
        ent[0] = eng[0];
        ent[1] = eng[1];
        ent[2] = eng[2];
        ent[3] = eng[2];
        ent[4] = eng[3];
        ent[5] = eng[4];
        ent[6] = eng[5];
        ent[7] = eng[5];
      }
    }
    else if (z.t==PYRAMID4) { //*3
      // these are represented by a FEBRICK, with coalesced nodes
      for (unsigned c=0; c<ve.size(); ++c) {
        const vector< unsigned >  eng = ve[c];  // element nodes, original (gambit)
              vector< unsigned >& ent = ve[c];  // ...,           modified (tecplot)
        ent.resize(8);
        ent[0] = eng[0];
        ent[1] = eng[1];
        ent[2] = eng[3];
        ent[3] = eng[2];
        ent[4] = eng[4];
        ent[5] = eng[4];
        ent[6] = eng[4];
        ent[7] = eng[4];
      }
    }
    else if (z.t==FEBRICK) { //*1
      // these need renumbering
      for (unsigned c=0; c<ve.size(); ++c) {
        swap(ve[c][2],ve[c][3]);
        swap(ve[c][6],ve[c][7]);
      }
    }
    const int sharefrom = (type==SmURF::ORDERED || !i? -1:0);
    mwriter.writeZoneData(type,pack,ve,m.vv,sharefrom);
  }
  if (!m.z() && m.v()) {
    // no connectivities present, but there is a point cloud
    mwriter.writeZoneData(SmURF::ORDERED,pack,vector< vector< unsigned > >(),m.vv,-1);
  }
}


const vector< vector< unsigned > > f_smurf::convert_from_vtelem(const vector< melem >& ve1)
{
  vector< vector< unsigned > > ve2(ve1.size());
  for (unsigned j=0; j<ve2.size(); ++j)
    ve2[j] = ve1[j].n;
  return ve2;
}

const vector< melem > f_smurf::convert_to_vtelem(const vector< vector< unsigned > >& ve1)
{
  vector< melem > ve2(ve1.size());
  for (unsigned j=0; j<ve2.size(); ++j)
    ve2[j].n = ve1[j];
  return ve2;
}

