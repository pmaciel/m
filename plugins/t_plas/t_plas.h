#ifndef t_plas_h
#define t_plas_h

#include "mkernel.h"
#include "plas.h"


// module to use the Particle Lagrangian Solver (PLaS)
class t_plas : public m::mtransform,
               public plas {
 public:
  void transform(GetPot& o, m::mmesh& m);

 private:  // plas interface functions
   void setFlowSolverParamOnInit(PLAS_FLOWSOLVER_PARAM *fp);
   void getElmNodes(const int iz, const int ie, int *enodes) {
     const std::vector< unsigned >& n = M->vz[iz].e2n[ie].n;
     for (size_t i=0; i<n.size(); ++i)
       enodes[i] = (int) n[i];
   }
   plas_elmtype_t getElmType(const int iz, const int ie) {
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
   double getEulerianTimeScale(int nod)                           { return 0.; }
   double getOldQuantity      (const plas_quantity_t& Q, int nod) { return getQuantity(Q,nod); }
   double getQuantity         (const plas_quantity_t& Q, int nod) { return (m_quantity_idx[Q]<0? 0. : ( M->vv[ m_quantity_idx[Q] ][nod] )); }
   int    getWallBndFlag      (int bnd)                           { return m_ziswall[bnd]? 1:0; }
   void   screenOutput        (const std::string& text)           { std::cout << "t_plas: info: " << text << std::endl; }
   void   screenWarning       (const std::string& text)           { std::cout << "t_plas: warn: " << text << std::endl; }

 private:  // member variables
  m::mmesh *M;
  std::vector< int  > m_quantity_idx;
  std::vector< bool > m_ziswall;
  double m_dt;
  int m_iter;
  int m_numiter;

};


#endif
