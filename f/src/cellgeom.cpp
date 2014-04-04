
/* computes scaled inward normals and volume of cell */

#include "common.h"

void cellgeom(int ic, local_node_struct *No_loc, double *vol, int *inc_min)
{
  int id;
  int inc;
  int j1;
  int j2;
  double ndotr;
  double sign;
  double rvec[2];
  double d03[3];
  double d13[3];
  double d23[3];
  double d01[3];
  double d21[3];
  double d31[3];

  const std::vector< unsigned >& inode(e2n[ic].n);
  for ( inc=0 ; inc<Nvtcell ; inc++ ) {
    No_loc[inc].node = inode[inc] ;
    for ( id=0 ; id<Ndim ; id++ )
      No_loc[inc].norm[id] = 0. ;
  }
  *inc_min=0 ;

  if ( Ndim==2 ) {
    /* triangles */

    for ( inc=0 ; inc<Nvtcell ; inc++ ) {
      j1 = (inc+1)%Nvtcell ;
      j2 = (inc+2)%Nvtcell ;

      /* components of scaled normal */
      No_loc[inc].norm[0] = M.vv[1][inode[j1]] - M.vv[1][inode[j2]] ;
      No_loc[inc].norm[1] = M.vv[0][inode[j2]] - M.vv[0][inode[j1]] ;

      /* components of vector inc->j2 */
      rvec[0] = M.vv[0][inode[inc]] - M.vv[0][inode[j2]] ;
      rvec[1] = M.vv[1][inode[inc]] - M.vv[1][inode[j2]] ;

      /* check normal is inward-facing */
      ndotr = No_loc[inc].norm[0]*rvec[0] + No_loc[inc].norm[1]*rvec[1] ;
      if ( ndotr<0. ) {
        No_loc[inc].norm[0] = -No_loc[inc].norm[0] ;
        No_loc[inc].norm[1] = -No_loc[inc].norm[1] ;
      } ;

      /* square of magnitude */
      No_loc[inc].norm2 = No_loc[inc].norm[0]*No_loc[inc].norm[0] + No_loc[inc].norm[1]*No_loc[inc].norm[1] ;
      if ( No_loc[inc].norm2<No_loc[*inc_min].norm2 )
        *inc_min=inc ;

    } ;

  }

  else {
    /* tetrahedra */

    /* vectors joining nodes of cell */
    for ( id=0 ; id<Ndim ; id++ ) {
      d23[id] = M.vv[id][inode[2]] - M.vv[id][inode[3]] ;
      d03[id] = M.vv[id][inode[0]] - M.vv[id][inode[3]] ;
      d13[id] = M.vv[id][inode[1]] - M.vv[id][inode[3]] ;

      d01[id] = M.vv[id][inode[0]] - M.vv[id][inode[1]] ;
      d21[id] = M.vv[id][inode[2]] - M.vv[id][inode[1]] ;
      d31[id] = -d13[id] ;
    } ;

    /* vector products (either all inward or all outward) */
    vecprd3(d23,d13,No_loc[0].norm) ;
    vecprd3(d03,d23,No_loc[1].norm) ;
    vecprd3(d01,d31,No_loc[2].norm) ;
    vecprd3(d21,d01,No_loc[3].norm) ;

    /* halve magnitude and set orientation to inward */
    for ( id=0, sign=0. ; id<Ndim ; id++ )
      sign += d03[id]*No_loc[0].norm[id] ;
    sign = 0.5*(sign>=0.? 1.:-1.);
    for ( inc=0 ; inc<Nvtcell ; inc++ )
      for ( id=0 ; id<Ndim ; id++ )
        No_loc[inc].norm[id] *= sign ;

    for ( inc=0 ; inc<Nvtcell ; inc++ ) {
      No_loc[inc].norm2 =
        No_loc[inc].norm[0]*No_loc[inc].norm[0] +
        No_loc[inc].norm[1]*No_loc[inc].norm[1] +
        No_loc[inc].norm[2]*No_loc[inc].norm[2];
      if ( No_loc[inc].norm2 < No_loc[*inc_min].norm2 )
        *inc_min=inc ;
    } ;

  } ;

  /* cell volume */
  for ( inc=0, *vol=0. ; inc<Nvtcell ; inc++ )
    for ( id=0 ; id<Ndim ; id++ )
      *vol += M.vv[id][inode[inc]]*No_loc[inc].norm[id] ;
  *vol = *vol/((double)Ndim*(double)Nvtfce) ;

}

