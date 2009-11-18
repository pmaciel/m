
#include "mfactory.h"
#include "t_math.h"

using namespace std;
using namespace m;


Register< mtransform,t_math > mt_math(10,"-tRx","[d] x axis rotation angle, in radians",
                                         "-tRy","[d] y axis ...",
                                         "-tRz","[d] z axis ...",
                                         "-tS", "[d] x,y,z axis scaling factor",
                                         "-tSx","[d] x axis scaling factor",
                                         "-tSy","[d] y axis ...",
                                         "-tSz","[d] z axis ...",
                                         "-tTx","[d] x axis translation",
                                         "-tTy","[d] y axis ...",
                                         "-tTz","[d] z axis ...");


void t_math::transform(GetPot& o, mmesh& m)
{
  // get dimensions, key and value
  const unsigned dim = m.d();
  const string k = o[o.get_cursor()];
  const double v = (k.find("-tR")==0? o.get(o.inc_cursor(),0.) :
                   (k.find("-tS")==0? o.get(o.inc_cursor(),1.) :
                   (k.find("-tT")==0? o.get(o.inc_cursor(),0.) : 0. )));

  // find vector fields, for rotation, creating:
  // x,y,z groups of
  //   lists of
  //     pointers to the values vector of each variable
  vector< vector< vector< double >* > > vvf(dim);
  if ((k=="-tRx" && dim>2) || (k=="-tRy" && dim>2) || (k=="-tRz" && dim>1)) {
    const vector< bool > visvector = m.vvectors();
    for (unsigned i=dim; i<visvector.size(); ++i) {
      if (visvector[i])
        for (unsigned j=0; j<dim; ++j)
          vvf[j].push_back(&m.vv[i+j]);
    }
  }

  // apply
  if      (k=="-tRx" && dim>2)  rotate(m.vv[1],m.vv[2],vvf[1],vvf[2],v);
  else if (k=="-tRy" && dim>2)  rotate(m.vv[2],m.vv[0],vvf[2],vvf[0],v);
  else if (k=="-tRz" && dim>1)  rotate(m.vv[0],m.vv[1],vvf[0],vvf[1],v);

  else if (k=="-tS") {
    if (dim>0)  scale(m.vv[0],v);
    if (dim>1)  scale(m.vv[1],v);
    if (dim>2)  scale(m.vv[2],v);
  }
  else if (k=="-tSx" && dim>0)  scale(m.vv[0],v);
  else if (k=="-tSy" && dim>1)  scale(m.vv[1],v);
  else if (k=="-tSz" && dim>2)  scale(m.vv[2],v);

  else if (k=="-tTx" && dim>0)  translate(m.vv[0],v);
  else if (k=="-tTy" && dim>1)  translate(m.vv[1],v);
  else if (k=="-tTz" && dim>2)  translate(m.vv[2],v);
}


void t_math::rotate(vector< double >& vx, vector< double >& vy, vector< vector< double >* >& vvx, vector< vector< double >* >& vvy, const double a)
{
  // coordinates
  const double x0 = (*max_element(vx.begin(),vx.end()) - *min_element(vx.begin(),vx.end()))/2.;
  const double y0 = (*max_element(vy.begin(),vy.end()) - *min_element(vy.begin(),vy.end()))/2.;
  translate(vx,-x0);
  translate(vy,-y0);
  for (unsigned n=0; n<vx.size(); ++n) {
    const double x = vx[n];
    const double y = vy[n];
    vx[n] = x*cos(a) + y*sin(a);
    vy[n] =-x*sin(a) + y*cos(a);
  }
  translate(vx,x0);
  translate(vy,y0);

  // states
  for (unsigned d=0; d<vvx.size(); ++d) {
    vector< double >& vx = *vvx[d];
    vector< double >& vy = *vvy[d];
    for (unsigned n=0; n<vx.size(); ++n) {
      const double x = vx[n];
      const double y = vy[n];
      vx[n] = x*cos(a) + y*sin(a);
      vy[n] =-x*sin(a) + y*cos(a);
    }
  }
}


void t_math::scale(vector< double >& v, const double f)
{
  for (unsigned n=0; n<v.size(); ++n)
    v[n] *= f;
}


void t_math::translate(vector< double >& v, const double d)
{
  for (unsigned n=0; n<v.size(); ++n)
    v[n] += d;
}
