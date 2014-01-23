#ifndef t_pillaz_h
#define t_pillaz_h

#include "mkernel.h"
#include "pillaz.h"


// module to use the Particle Lagrangian Solver (Pillaz)
class t_pillaz : public m::mtransform,
               public pillaz {
 public:
  void transform(GetPot& o, m::mmesh& m, const XMLNode& x, const XMLNode& x);

 private:  // pillaz interface functions
   void setFlowSolverParamOnInit(PILLAZ_FLOWSOLVER_PARAM *_fp);
   void getBndElmNodes(const int iz, const int ie, int *enodes) { getElmNodes(iz,ie,enodes); }
   pillaz_elmtype_t getBndElmType(const int iz, const int ie) { return getElmType(iz,ie); }
   void getElmNodes(const int iz, const int ie, int *enodes) {
     const std::vector< unsigned >& n = M->vz[iz].e2n[ie].n;
     for (size_t i=0; i<n.size(); ++i)
       enodes[i] = (int) n[i];
   }
   pillaz_elmtype_t getElmType(const int iz, const int ie) {
     const m::mtype t(M->vz[iz].t);
     return (t==m::FELINESEG?       ELM_EDGE        :
            (t==m::FETRIANGLE?      ELM_TRIANGLE    :
            (t==m::FEQUADRILATERAL? ELM_QUAD        :
            (t==m::FETETRAHEDRON?   ELM_TETRAHEDRON :
            (t==m::FEBRICK?         ELM_BRICK       :
            (t==m::PRISM3?          ELM_WEDGE       :
            (t==m::PYRAMID4?        ELM_PYRAMID     :
                                    ELM_UNDEFINED )))))));
   }
   double getEulerianTimeScale(int nod){ return 0.; }
   double getQuantityScalar(const pillaz_quantityscalar_t& Q, const int nod) {
     return (m_quantityscalar[Q]<0? 0. : M->vv[m_quantityscalar[Q]][nod]);
   }
   void getQuantityVector(const pillaz_quantityvector_t& Q, const int nod, const int dim, double *v) {
     for (int i=0; i<dim; ++i)
       v[i] = (m_quantityvector[Q]<0? 0. : M->vv[m_quantityvector[Q]+i][nod]);
   }
   int    getWallBndFlag(int bnd)               { return m_ziswall[bnd]? 1:0; }
   void   screenOutput(const std::string& text) { std::cout << "t_pillaz: info: " << text << std::endl; }

 private:  // member variables
  m::mmesh *M;
  std::vector< int  > m_quantityscalar;
  std::vector< int  > m_quantityvector;
  std::vector< bool > m_ziswall;
  double m_dt;
  int m_numiter;

};


#endif
