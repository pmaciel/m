
/* read grid and starting solution from tecplot file */

#include <algorithm>
#include <memory>
#include "boost/progress.hpp"
#include "common.h"


/*
 * comparison of boundary element with "inner cells" connectivity, returning the
 * "inner cell" index and respective opposite node global/local indices.
 * in: "inner" connectivity
 * in: boundary element to compare
 * out: "inner" correponding cell (0-based)
 * out: ... cell opposite node, global index (0-based)
 * out: ... cell opposite node, local index (0-based)
 */
void fab(
  const std::vector< m::melem >& e2n, const std::vector< unsigned >& ben,
  int *bcell, int *bnode, int *binc )
{
  bool match = false;
  for (unsigned c=0; c<e2n.size() && !match; ++c) {
    const std::vector< unsigned >& cell = e2n[c].n;

    // check for "inner" element matching all nodes from the boundary element
    int n_match = 0;
    for (unsigned i=0; i<ben.size(); ++i)
      n_match += std::count(cell.begin(),cell.end(),ben[i])? 1:0;
    match = n_match==Nvtfce;

    // if found, locate which node index (local/global) is not contained
    if (match) {
      *bcell = (int) c;
      for (unsigned i=0; i<cell.size(); ++i)
        if (!std::count(ben.begin(),ben.end(),cell[i])) {
          *bnode = (int) cell[i];
          *binc  = (int)      i;
          break;
        }
    }
  }
  if (!match)
    nrerror("boundary element not found");
}


/*
 * read grid/solution file in tecplot format
 * (first into m:mmesh format, then convert into native data structure)
 * in: file name
 * in: if solution is to be read as well
 */
void readsoltp(const std::string& infile, int read_soln)
{
  using std::cout;
  using std::endl;
  using std::vector;


  cout << "readsoltp: reading m::mmesh from \"" << infile << "\"..." << endl;
  {
    std::auto_ptr< m::mfinput > p(m::mfactory< m::mfinput  >::instance()->Create(extension(infile)));
    char* argv[] = { (char*) "", const_cast< char* >(infile.c_str()) };
    GetPot o2(2,argv);
    p->read(o2,M);

    cout << "readsoltp: m::mmesh:  d=" << M.d()
              << "  n=" << M.n()
              << "  v=" << M.v()
              << endl;
    for (unsigned t=0; t<M.vz.size(); ++t)
      cout << "  zone:  d=" << M.d(t)
                << "  e=" << M.e(t)
                << "  n=" << M.vz[t].n
                << "  t=" << M.vz[t].t
                << endl;
    cout << "  variables: ";
    for (unsigned i=0; i<M.v(); ++i)
      cout << " \"" << M.vn[i] << '"';
    cout << endl;
  }

  // move out "inner" connectivity
  e2n.swap(M.vz[0].e2n);

  // move out wall distance
  No_wd.clear();
  if (M.vn.back()=="wd") {
    if (walldist)
      No_wd.swap(M.vv.back());  // (always the last variable)
    M.vv.pop_back();
    M.vn.pop_back();
  }
  cout << "readsoltp: reading m::mmesh." << endl;


  cout << "readsoltp: set mesh sizes..." << endl;
  Nnode = (int) M.n();
  Ncell = (int) e2n.size();
  const int Nvar = (int) M.v() - (int) M.d();
  if (read_soln && Nvar!=Nsys)
    nrerror("Not enough variables in solution file for restart");
  if (!read_soln) {
    M.vn.resize(Ndim+Nsys);
    M.vv.resize(Ndim+Nsys);
    for (int i=0; i<Nsys; ++i) {
      M.vn[Ndim+i] = m_vars_label[i];
      M.vv[Ndim+i].assign(Nnode,m_vars_init[i]);
    }
  }
  cout << "readsoltp: set mesh sizes." << endl;


  cout << "readsoltp: copy/initialize solution..." << endl;
  No_W.assign(Nsys,std::vector< double > (Nnode,0.));
  for (int i=0; i<Nsys; ++i)
    for (int n=0; n<Nnode; ++n)
      No_W[i][n] = M.vv[Ndim+i][n];
  cout << "readsoltp: copy/initialize solution." << endl;


  // count boundary elements
  Nbface = 0;
  for (unsigned i=1; i<M.z(); ++i)
    Nbface += (int) M.e(i);


  if (Nbface) {
    // allocate structures
    Fab_cell .resize(Nbface);
    Fab_node .resize(Nbface);
    Fab_inc  .resize(Nbface);
    Fab_group.resize(Nbface);

    // check for the boundary-node connectivity file and generate or read it
    const std::string infilefab = infile + ".fab";
    std::fstream f(infilefab.c_str(),std::ios::in);
    if (!f) {
      cout << "readsoltp: generate boundary-node connectivity (\"" << infilefab << "\")..." << endl;
      f.open(infilefab.c_str(),std::ios::out|std::ios_base::trunc);

      boost::progress_display pbar(Nbface);
      f << Nbface << endl;
      int ifc = 0;
      for (unsigned i=1; i<M.z(); ++i)
        for (unsigned j=0; j<M.e(i); ++j, ++ifc, ++pbar) {
          fab(
            e2n, M.vz[i].e2n[j].n,
            &Fab_cell[ifc], &Fab_node[ifc], &Fab_inc[ifc]);
          Fab_group[ifc] = (int) i;
          f << Fab_group[ifc] << ' '
            << Fab_cell [ifc] << ' '
            << Fab_node [ifc] << ' '
            << Fab_inc  [ifc] << endl;
        }

      f.close();
      cout << "readsoltp: generate boundary-node connectivity." << endl;
    }
    else {
      cout << "readsoltp: read boundary-node connectivity (\"" << infilefab << "\")..." << endl;
      int ifc = 0;
      f >> ifc;
      if (ifc!=Nbface)
        nrerror("The boundary-node file does not correspond to this mesh");
      boost::progress_display pbar(Nbface);
      for (ifc=0; ifc<Nbface && f; ++ifc, ++pbar)
        f >> Fab_group[ifc] >> Fab_cell[ifc] >> Fab_node[ifc] >> Fab_inc[ifc];
      if (!f)
        nrerror("The boundary-node file is badly-formed");
      cout << "readsoltp: read boundary-node connectivity." << endl;
    }

  }


  if (periodic) {
    cout << "readsoltp: construct periodic arrays..." << endl;

    // mark "periodic" nodes
    vector< int > No_peri(Nnode,-1);
    vector< int > No_twin(Nnode,-1);
    for (int f=0; f<Nbface; ++f) {
      const int IB = BCgroup[Fab_group[f]].type;
      if (IB==IBPERI || IB==IBPERE) {
        const int c = Fab_cell[f];
        for (int i=0; i<Nvtcell; ++i)
          if ((int) (e2n[c].n[i])!=Fab_node[f])
            No_peri[ e2n[c].n[i] ] = IB;
      }
    }
    const int Nper = Nnode - (int) std::count(No_peri.begin(),No_peri.end(),-1);
    if (Nper%2)
      nrerror("odd number of periodic nodes");

    // set twin nodes
    for (int n=0; n<Nnode; ++n) {
      if (No_peri[n]==IBPERI) {
        for (int m=0; m<Nnode; ++m) {
          if (No_peri[m]==IBPERE) {

            int flag = 0;
            for (int id=1; id<Ndim; id++) {
              const int jd = (id+periodic_dirn)%Ndim;
              if (std::abs(M.vv[jd][n] - M.vv[jd][m])<1.e-10)
                flag++;
            }
            if (flag==Ndim-1) {
              No_twin[n] = m;
              No_twin[m] = n;
              break;
            }
          }

        }
        if (No_twin[n]<0)
          nrerror("No twin found for periodic node");
      }
    }

    e2n_periodic = e2n;  // copy!
    for (unsigned i=0; i<(unsigned) e2n.size(); ++i) {
      const vector< unsigned >& ne_normal   = e2n[i].n;
            vector< unsigned >& ne_periodic = e2n_periodic[i].n;
      for (unsigned j=0; j<(unsigned) ne_normal.size(); ++j)
        if (No_peri[ ne_normal[j] ]==IBPERE)
          ne_periodic[j] = No_twin[ ne_normal[j] ];
    }

    cout << "readsoltp: construct periodic arrays." << endl;
  }


  cout << "readsoltp: setup sparsity..." << endl;
  vector< vector< unsigned > > nz(Nnode);
  for (unsigned c=0; c<(unsigned) e2n.size(); ++c) {
    const vector< unsigned >& en = e2n[c].n;
    for (unsigned i=0; i<(unsigned) en.size(); ++i)
      for (unsigned j=0; j<(unsigned) en.size(); ++j) {
        nz[ en[i] ].push_back( en[j] );
        nz[ en[j] ].push_back( en[i] );
        if (periodic) {
          const vector< unsigned >& en = e2n_periodic[c].n;
          nz[ en[i] ].push_back( en[j] );
          nz[ en[j] ].push_back( en[i] );
        }
      }
  }
  for (vector< vector< unsigned > >::iterator n=nz.begin(); n!=nz.end(); ++n) {
    std::sort(n->begin(),n->end());
    n->erase(std::unique(n->begin(),n->end()),n->end());
  }
  cout << "readsoltp: setup sparsity." << endl;


  if (temperature && !scalar_coupling) {
    cout << "readsoltp: ls_scalar setup..." << endl;
    ls_scalar->initialize(Nnode,Nnode,1);
    if (ls_scalar->issparse) ls_scalar->initialize(nz);
    cout << "readsoltp: ls_scalar setup." << endl;
  }
  if (turmod && turbulence_coupling==0) {
    cout << "readsoltp: ls_turb{1,2} setup..." << endl;
    ls_turb1->initialize(Nnode,Nnode,1);
    ls_turb2->initialize(Nnode,Nnode,1);
    if (ls_turb1->issparse) ls_turb1->initialize(nz);
    if (ls_turb2->issparse) ls_turb2->initialize(nz);
    cout << "readsoltp: ls_turb{1,2} setup." << endl;
  }
  if (turmod && turbulence_coupling==1) {
    cout << "readsoltp: ls_turb setup..." << endl;
    ls_turb->initialize(Nnode,Nnode,2);
    if (ls_turb->issparse) ls_turb->initialize(nz);
    cout << "readsoltp: ls_turb setup." << endl;
  }
  if (true) {
    cout << "readsoltp: ls_coupled setup..." << endl;
    ls_coupled->initialize(Nnode,Nnode,Ncoupled);
    if (ls_coupled->issparse) ls_coupled->initialize(nz);
    cout << "readsoltp: ls_coupled setup." << endl;
  }


  cout << "readsoltp: assigning boundary-node structure(s)..." << endl;

  // set boundary-node structure
  Nobg = new bound_node_struct*[1+Nbcgroup];
  for (int g=1; g<=Nbcgroup; ++g) {

    // mark which nodes sit in this boundary
    vector< int > gnodes;
    gnodes.reserve(Nnode);
    for (int f=0; f<Nbface; ++f) {
      const int c = Fab_cell[f];
      for (int i=0; i<Nvtcell && Fab_group[f]==g; ++i)
        if (i!=Fab_inc[f])
          gnodes.push_back((int) e2n[c].n[i]);
    }

    // count nodes and resize Nobg (and switch off groups with no faces)
    BCgroup[g].nnode = (int) gnodes.size();
    BCgroup[g].type = gnodes.size()? BCgroup[g].type : IBNONE;
    if (!BCgroup[g].type)
      continue;
    Nobg[g] = new bound_node_struct[1+BCgroup[g].nnode];

    // set Nobg nodes
    for (int i=0; i<(int) gnodes.size(); ++i) {
      Nobg[g][1+i].node = gnodes[i];
      Nobg[g][1+i].modn = 0.;
      for (int d=0; d<Ndim; ++d)
        Nobg[g][1+i].n[d] = 0.;
    }

    // find normal for each boundary group
    int bf = -1;
    for (int f=0; f<Nbface && bf<0; ++f)
      bf = Fab_group[f]==g? f:bf;
    if (bf<0)
      nrerror("boundary face for this group wasn't found");

    local_node_struct No_local[4];
    int inc_min;
    double vol;
    cellgeom(Fab_cell[bf],No_local,&vol,&inc_min);

    // normalize it and assign
    const int o = Fab_inc[bf];   // local opposite node index
    double modn = 0.;
    for (int d=0; d<Ndim; ++d)
      modn += No_local[o].norm[d]*No_local[o].norm[d];
    modn = sqrt(modn);
    for (int d=0; d<Ndim; ++d)
      BCgroup[g].n[d] = No_local[o].norm[d]/modn;

  }

  // mark nodes boundary groups, assigned according to the boundary index
  // (with priority to walls)
  No_group.assign(Nnode,0);
  for (int g=1; g<=Nbcgroup; ++g) {
    if (BCgroup[g].type)
      for (int inb=1; inb<=BCgroup[g].nnode; ++inb)
        No_group[Nobg[g][inb].node] = g;
  }
  for (int g=1; g<=Nbcgroup; ++g) {
    if (BCgroup[g].type==IBWALL || BCgroup[g].type==IBWALQ)
      for (int inb=1; inb<=BCgroup[g].nnode; ++inb)
        No_group[Nobg[g][inb].node] = g;
  }

  cout << "readsoltp: assigning boundary-node structure(s)." << endl;


  cout << "readsoltp: allocate memory for nodal structures..." << endl;
  No_vol.assign(Nnode,0.);
  if (turmod) {
    No_nuturb.assign(Nnode,0.);
    No_dissipation.assign(Nnode,0.);
    if (walldist)
      No_lenturb.assign(Nnode,0.);
  }
  cout << "readsoltp: allocate memory for nodal structures." << endl;


  if (turmod && !read_soln) {
    cout << "readsoltp: set initial turbulence values..." << endl;
    const double turb_intensity = m_vars_init[iv_turb1];
    const double turb_width     = m_vars_init[iv_turb2];  // (reference length)

    // set initial k
    double q2 = 0.;
    for (int d=1; d<=Ndim; ++d)
      q2 += m_vars_init[d]*m_vars_init[d];
    const double k = 1.5*turb_intensity*turb_intensity*q2;

    // set initial epsilon/omega
    const double e = turmod/10==ITMGKE?
      pow(Cmu,0.75)*pow(k,1.5)/(0.09*turb_width) : // epsilon
      sqrt(k)/(0.09*turb_width);                   // omega

    for (int n=0; n<Nnode; ++n) {
      No_W[iv_turb1][n] = k;
      No_W[iv_turb2][n] = e;
    }

    // initialize turbulence length
    if (walldist) {
      // ... in the field
      for (int n=0; n<Nnode; ++n) {
        const double len = pow(1. - No_wd[n]/turb_width,2.);
        No_lenturb[n] = 0.53*turb_width*(0.14-len*(0.08+0.06*len));
      }
      // ... at the walls (close to zero)
      for (int f=0; f<Nbface; ++f) {
        if (BCgroup[Fab_group[f]].type==IBWALL || BCgroup[Fab_group[f]].type==IBWALQ)
          for (int inc=0; inc<Nvtcell; ++inc)
            if (inc!=Fab_inc[f])
              No_lenturb[ e2n[Fab_cell[f]].n[inc] ] = 1.e-20;
      }
    }

    cout << "readsoltp: set initial turbulence values." << endl;
  }


  cout << "initializing boundary conditions: most of them..." << endl;
  for (int ig=1; ig<=Nbcgroup; ig++) {

  int invals_offset = 0;  // for IBFIXV only
  switch (BCgroup[ig].type) {
    case IBNONE:                                 /* GROUP NOT ACTIVE */
      break;

    case IBFIXV:                                 /* FIXED FLUX */

      /* Uniform velocity */
      if (BCgroup[ig].option==0) {
        invals_offset=0;
        for (int inb=1; inb<=BCgroup[ig].nnode; inb++) {
          const int inu = Nobg[ig][inb].node;
          for (int iv=1; iv<=Ndim; iv++)
            No_W[iv][inu]=BCgroup[ig].invals[0]*BCgroup[ig].n[iv-1];
        }
      }

      /* Parabolic velocity profile */
      else if (BCgroup[ig].option==1) {
        invals_offset=2;

        const int idirn = std::abs(BCgroup[ig].n[0])>=std::abs(BCgroup[ig].n[1])? 1:0;

        double xymax=-1.e10, xymin=1.e10;
        int    inmin=0,      inmax=0;
        for (int inb=1; inb<=BCgroup[ig].nnode; inb++) {
          if (M.vv[idirn][Nobg[ig][inb].node]<xymin) {
            inmin=Nobg[ig][inb].node;
            xymin=M.vv[idirn][inmin];
          }
          if (M.vv[idirn][Nobg[ig][inb].node]>xymax) {
            inmax=Nobg[ig][inb].node;
            xymax=M.vv[idirn][inmax];
          }
        }

        double dinlet=0.;
        for (int id=0; id<Ndim; id++)
          dinlet+=(M.vv[id][inmax]-M.vv[id][inmin])*(M.vv[id][inmax]-M.vv[id][inmin]);
        dinlet = sqrt(dinlet);
        const double d = (BCgroup[ig].invals[2]-BCgroup[ig].invals[1])*dinlet;

        /* calculate maximum velocity from volume-averaged value */
        const double Vmax=1.5*BCgroup[ig].invals[0];

        for (int inb=1; inb<=BCgroup[ig].nnode; inb++) {
          const int inu = Nobg[ig][inb].node;
          double r=0.;
          for (int id=0; id<Ndim; id++)
            r += (M.vv[id][inu]-M.vv[id][inmin])*(M.vv[id][inu]-M.vv[id][inmin]);
          r = (sqrt(r) - BCgroup[ig].invals[1]*dinlet)/d;
          for (int iv=1; iv<=Ndim; iv++)
            No_W[iv][inu] = 4.0*Vmax*r*(1.-r)*BCgroup[ig].n[iv-1];
        }
      }

      /* Parabolic velocity profile in circular duct */
      else if (BCgroup[ig].option==4) {
        invals_offset=4;
        std::vector< double > Xcen(Ndim,0.);
        const double Vmax=2.0*BCgroup[ig].invals[0];

        for (int id=0; id<Ndim; id++)
          Xcen[id]=BCgroup[ig].invals[id+1];

        const double Radius = BCgroup[ig].invals[Ndim+1];

        for (int inb=1; inb<=BCgroup[ig].nnode; inb++) {
          const int inu = Nobg[ig][inb].node;
          double r=0.;
          for (int id=0; id<Ndim; id++)
            r += (M.vv[id][inu]-Xcen[id])*(M.vv[id][inu]-Xcen[id]);
          r = sqrt(r)/Radius;
          for (int iv=1; iv<=Ndim; iv++)
            No_W[iv][inu] = Vmax*(1.-r*r)*BCgroup[ig].n[iv-1];
        }
      }

      /* Boundary-layer velocity profile */
      else if (BCgroup[ig].option==2) {
        invals_offset=2;

        const int idirn = std::abs(BCgroup[ig].n[0])>=std::abs(BCgroup[ig].n[1])? 1:0;

        double xymax=-1.e10, xymin=1.e10;
        int    inmin=0,      inmax=0;
        for (int inb=1; inb<=BCgroup[ig].nnode; inb++) {
          if (M.vv[idirn][Nobg[ig][inb].node]<xymin) {
            inmin=Nobg[ig][inb].node;
            xymin=M.vv[idirn][inmin];
          }
          if (M.vv[idirn][Nobg[ig][inb].node]>xymax) {
            inmax=Nobg[ig][inb].node;
            xymax=M.vv[idirn][inmax];
          }
        }

        double dinlet=0.;
        for (int id=0; id<Ndim; id++)
          dinlet+=(M.vv[id][inmax]-M.vv[id][inmin])*(M.vv[id][inmax]-M.vv[id][inmin]);
        dinlet = sqrt(dinlet);

        const double Ue=BCgroup[ig].invals[0];
        const double deltar=BCgroup[ig].invals[1]*dinlet;
        const double deltal=BCgroup[ig].invals[2]*dinlet;

        double x = pow(Ue/nulam,0.25)*pow(2.70*deltar,1.25);
        double Rx = Ue*x/nulam;
        double utaur = sqrt(0.5*0.0592*Ue*Ue*pow(Rx,-0.2));
        const double Cr = (Ue/utaur)*pow(utaur*deltar/nulam,-0.143);

        x = pow(Ue/nulam,0.25)*pow(2.70*deltal,1.25);
        Rx = Ue*x/nulam;
        double utaul = sqrt(0.5*0.0592*Ue*Ue*pow(Rx,-0.2));
        const double Cl = (Ue/utaul)*pow(utaul*deltal/nulam,-0.143);

        for (int inb=1; inb<=BCgroup[ig].nnode; inb++) {
          const int inu = Nobg[ig][inb].node;
          double y=0.;
          for (int id=0; id<Ndim; id++)
            y += (M.vv[id][inu]-M.vv[id][inmin])*(M.vv[id][inu]-M.vv[id][inmin]);
          y = sqrt(y);
          double Uy = -1.;
          if (y<deltar) {
            if ((y*utaur/nulam)>11.8)
              Uy = Cr*utaur*pow(utaur*y/nulam,0.143);
            else
              Uy = utaur*utaur*y/nulam;
          }
          else if (y>(dinlet-deltal)) {
            if (((dinlet-y)*utaul/nulam)>11.8)
              Uy = Cl*utaul*pow(utaul*(dinlet-y)/nulam,0.143);
            else
              Uy = utaul*utaul*(dinlet-y)/nulam;
          }
          else
            Uy = Ue;
          for (int iv=1; iv<=Ndim; iv++)
            No_W[iv][inu] = Uy*BCgroup[ig].n[iv-1];
        }
      }

      /* Read inlet profiles from inlet file */
      else if (BCgroup[ig].option==3) {
        if (Ndim==2)
          read_inlet_2D(file_inlet,ig);
        else if (Ndim==3)
          read_inlet_3D(file_inlet,ig);
      }

/* Set scalar and turbulence inlet values */

      if (temperature && BCgroup[ig].option!=3) {
        for (int inb=1; inb<=BCgroup[ig].nnode; inb++) {
          const int inu = Nobg[ig][inb].node;
          No_W[iv_temp][inu]=BCgroup[ig].invals[invals_offset+1];
        }
        invals_offset++;
      }

      if (turmod && BCgroup[ig].option!=3)
      for (int inb=1; inb<=BCgroup[ig].nnode; inb++) {
        const int inu = Nobg[ig][inb].node;

        const double turb_intensity = BCgroup[ig].invals[invals_offset+1];
        const double turb_width     = BCgroup[ig].invals[invals_offset+2];

        double q2=0.;
        for (int iv=1; iv<=Ndim; iv++)
          q2 += No_W[iv][inu]*No_W[iv][inu];
        const double k = turb_intensity*turb_intensity*q2;

        const double e = turmod/10==ITMGKE?   // epsilon (or omega)
                         Cmu*pow(k,1.5)/(0.07*turb_width)
           // pow(Cmu,3./4.)*pow(k,1.5)/(0.07*turb_width)
           //            Cmu*pow(k,2  )/(0.01*nulam)
                       : pow(k,0.5)/(0.07*turb_width);

        No_W[iv_turb1][inu] = k;
        No_W[iv_turb2][inu] = e;
      }

      break;

    case IBFIXP:                                 /* FIXED PRESSURE */

      for(int inb=1; inb<=BCgroup[ig].nnode; inb++)
        No_W[0][Nobg[ig][inb].node]=BCgroup[ig].invals[0];

      if (BCgroup[ig].option!=0) {
        for (int iv=1; iv<=Ndim; iv++)
          if (iv!=BCgroup[ig].option)
            for(int inb=1; inb<=BCgroup[ig].nnode; inb++)
              No_W[iv][Nobg[ig][inb].node]=0.;
      }

      break;

    case IBSYMM:                                 /* SYMMETRY PLANE (X, Y or Z DIRN.) */

      for(int inb=1; inb<=BCgroup[ig].nnode; inb++)
        No_W[ BCgroup[ig].option ][Nobg[ig][inb].node]=0.;

      break;

    default:
      break;

  }
  }
  cout << "initializing boundary conditions: most of them." << endl;


  cout << "initializing boundary conditions: moving walls..." << endl;
  for (int ig=1; ig<=Nbcgroup; ig++) {
    if ((BCgroup[ig].type==IBWALL || BCgroup[ig].type==IBWALQ) && BCgroup[ig].option==1) {
      /* wall speed components (u,v,w) specified */
      for (int inb=1; inb<=BCgroup[ig].nnode; inb++)
        for (int iv=1; iv<=Ndim; iv++)
          No_W[iv][Nobg[ig][inb].node]=BCgroup[ig].invals[iv-1+temperature];
    }
  }
  cout << "initializing boundary conditions: moving walls." << endl;


  cout << "initializing boundary conditions: stationary walls..." << endl;
  for (int ig=1; ig<=Nbcgroup; ig++)
    if ((BCgroup[ig].type==IBWALL || BCgroup[ig].type==IBWALQ) && BCgroup[ig].option==0)
    for (int inb=1; inb<=BCgroup[ig].nnode; inb++) {
      for (int iv=1; iv<=Ndim; iv++)
        No_W[iv][Nobg[ig][inb].node]=0.;
    if (turmod)
      No_W[iv_turb1][Nobg[ig][inb].node]=1.e-20;
  }
  cout << "initializing boundary conditions: stationary walls." << endl;


  if (temperature) {
    cout << "initializing boundary conditions: fixed-temperature walls..." << endl;
    for (int ig=1; ig<=Nbcgroup; ig++)
      if (BCgroup[ig].type==IBWALL)
        for (int inb=1; inb<=BCgroup[ig].nnode; inb++)
          No_W[iv_temp][Nobg[ig][inb].node]=BCgroup[ig].invals[0];
    cout << "initializing boundary conditions: fixed-temperature walls." << endl;
  }


  if (turmod) {
    cout << "initializing boundary conditions: turbulent k=0 at the walls..." << endl;
    for (int ig=1; ig<=Nbcgroup; ig++)
      if (BCgroup[ig].type==IBWALL || BCgroup[ig].type==IBWALQ) {
        for (int inb=1; inb<=BCgroup[ig].nnode; inb++)
          No_W[iv_turb1][Nobg[ig][inb].node]=1.e-20;
      }
    cout << "initializing boundary conditions: turbulent k=0 at the walls." << endl;
  }


  cout << "calculating cell geometry..." << endl;
  // carries out various geometric and counting operations on grid:
  //  boundary nodes, nodal (dual) volumes, domain volume, degree
  int inc_min;
  double vol;
  local_node_struct No_local[4];


/* Mark wall cells (containing at least one wall node) */
  Ce_type.assign(Ncell,0);
  for (int ic=0; ic<Ncell; ic++)
    for (int inc=0; inc<Nvtcell && !Ce_type[ic]; inc++) {
      const int IB = BCgroup[No_group[ e2n[ic].n[inc] ]].type;
      Ce_type[ic] = (IB==IBWALL || IB==IBWALQ)? 1:0;
    }


/* Calculate domain volume and dual volume and degree of nodes */
  cout << "calculating cell geometry: nodal volumes..." << endl;
/* Nodal volumes and grid quality */
/* N.B. These are further modified for periodic nodes at end of routine */
  for (int ic=0; ic<Ncell; ic++) {
    cellgeom(ic,No_local,&vol,&inc_min);

    for (int inc=0; inc<Nvtcell; inc++)
      No_vol[No_local[inc].node] += vol/(double)Nvtcell;
  }

  Vtot = 0.;
  for (int inu=0; inu<Nnode; inu++)
    Vtot += No_vol[inu];
  cout << "calculating cell geometry: volume=" << Vtot << endl;


  cout << "find nearest non-wall node (assumed to be neighbour)..." << endl;
  if (turmod && !wall_functions) {
    for (int ig=1; ig<=Nbcgroup; ig++)
      if (BCgroup[ig].type==IBWALL || BCgroup[ig].type==IBWALQ)
        for (int inb=1; inb<=BCgroup[ig].nnode; inb++) {

      const int inu = Nobg[ig][inb].node;
      Nobg[ig][inb].dist = 1.e20;
      Nobg[ig][inb].twin = inu;

      for (int i=0; i<(int) nz[inu].size(); ++i) {
        if ((int) nz[inu][i]==inu)
          continue;
        const int inhbr = (int) nz[inu][i];
        const int IB = BCgroup[No_group[inhbr]].type;
        if (IB!=IBWALL && IB!=IBWALQ) {
          double wdist = 0.;
          for (int id=0; id<Ndim; id++)
            wdist += (M.vv[id][inu]-M.vv[id][inhbr])*(M.vv[id][inu]-M.vv[id][inhbr]);
          wdist = sqrt(wdist);
          if (wdist<Nobg[ig][inb].dist) {
            Nobg[ig][inb].dist = wdist;
            Nobg[ig][inb].twin = inhbr;
          }
        }
      }

      /*
       * If all the neighbour nodes are wall nodes then search whole grid for
       * nearest non-wall node for twin (which isn't a neighbour)
       */
      if (Nobg[ig][inb].twin==inu) {
        Nobg[ig][inb].dist=1.e20;
        for (int jnu=0; jnu<Nnode; jnu++) {
          double wdist = 0.;
          for (int id=0; id<Ndim; id++)
            wdist += (M.vv[id][inu]-M.vv[id][jnu])*(M.vv[id][inu]-M.vv[id][jnu]);
          wdist = sqrt(wdist);
          const int IB = BCgroup[No_group[jnu]].type;
          if (IB!=IBWALL && IB!=IBWALQ)
          if (wdist<Nobg[ig][inb].dist) {
            Nobg[ig][inb].dist = wdist;
            Nobg[ig][inb].twin = jnu;
          }
        }
      }

        }
  }
  cout << "find nearest non-wall node (assumed to be neighbour)." << endl;


  cout << "calculate wall-node normals..." << endl;
  std::vector< std::vector< double > > No_norm(Ndim,std::vector< double >(Nnode,0.));

  for (int ifc=0; ifc<Nbface; ifc++)
  if (Ce_type[Fab_cell[ifc]]==1) {
    const int ic = Fab_cell[ifc];
    const int inc_op = Fab_inc[ifc];

    cellgeom(ic,No_local,&vol,&inc_min);

    for (int inf=1; inf<=Nvtfce; inf++) {
      const int inc=(inf+inc_op)%Nvtcell;
      const int inu=e2n[ic].n[inc];
      for (int id=0; id<Ndim; id++)
        No_norm[id][inu] += No_local[inc_op].norm[id]/dNvtfce;
    }
  }

  for (int ig=1; ig<=Nbcgroup; ig++)
    if (BCgroup[ig].type==IBWALL || BCgroup[ig].type==IBWALQ)
    for (int inb=1; inb<=BCgroup[ig].nnode; inb++) {
      const int inu = Nobg[ig][inb].node;
    Nobg[ig][inb].modn=0.;
    for (int id=0; id<Ndim; id++) {
      Nobg[ig][inb].n[id] = No_norm[id][inu];
      Nobg[ig][inb].modn += No_norm[id][inu]*No_norm[id][inu];
    }
    Nobg[ig][inb].modn=sqrt(Nobg[ig][inb].modn);
  }
  cout << "calculate wall-node normals." << endl;


  if (turmod && wall_functions) {
    cout << "wall functions initialization..." << endl;
    // WF nodes are more appropriate for high-Re turbulence models
    // (group number Nbcgroup+1)

    cout << "wall functions: mark and count wall-function nodes..." << endl;
    int Nwfnode = 0;
    for (int ic=0; ic<Ncell; ic++) {
      if (Ce_type[ic]==1)
        for (int inc=0; inc<Nvtcell; inc++) {
          const int inu = e2n[ic].n[inc];
          if ((No_group[inu]==0 || BCgroup[No_group[inu]].type==IBFIXP || BCgroup[No_group[inu]].type==IBPERI)) {
            No_group[inu]=Nbcgroup+1;
            Nwfnode++;
          }
        }
    }

    BCgroup[Nbcgroup+1].nnode=Nwfnode;
    cout << "wall functions: mark and count wall-function nodes." << endl;


    cout << "wall functions: find nearest wall node (assumed to be neighbour)..." << endl;
    WFnodes.resize(Nwfnode);  // allocate
    int inwf=-1;
    for (int inu=0; inu<Nnode; inu++)
      if (BCgroup[No_group[inu]].type==IBWLFN) {
        inwf++;
        WFnodes[inwf].node = inu;
        WFnodes[inwf].dist = 1.e20;

        for (int i=0; i<(int) nz[inu].size(); ++i) {
          if ((int) nz[inu][i]==inu)
            continue;
          const int inhbr = (int) nz[inu][i];
          const int IB = BCgroup[No_group[inhbr]].type;
          if (IB==IBWALL || IB==IBWALQ) {
            double wdist = 0.;
            for (int id=0; id<Ndim; id++)
              wdist += (M.vv[id][inu]-M.vv[id][inhbr])*(M.vv[id][inu]-M.vv[id][inhbr]);
            wdist = sqrt(wdist);
            if (wdist<WFnodes[inwf].dist) {
              WFnodes[inwf].dist = wdist;
              WFnodes[inwf].twin = inhbr;
            }
          }
        }
      }
    cout << "wall functions: find nearest wall node (assumed to be neighbour)." << endl;


    cout << "wall functions: calculate nodal normals at wall nodes..." << endl;
    No_norm.assign(Ndim,std::vector< double >(Nnode,0.));  // allocate

    for (int ifc=0; ifc<Nbface; ifc++)
    if (BCgroup[Fab_group[ifc]].type==IBWALL || BCgroup[Fab_group[ifc]].type==IBWALQ) {
      const int ic = Fab_cell[ifc];
      const int inc_op = Fab_inc[ifc];

      cellgeom(ic,No_local,&vol,&inc_min);

      // switch to periodic connectivity
      if (periodic)
        e2n.swap(e2n_periodic);

      for (int inf=1; inf<=Nvtfce; inf++) {
        const int inc=(inf+inc_op)%Nvtcell;
        const int inu=e2n[ic].n[inc];
        for (int id=0; id<Ndim; id++)
          No_norm[id][inu] += No_local[inc_op].norm[id]/dNvtfce;
      }

      // switch to "normal" connectivity
      if (periodic)
        e2n.swap(e2n_periodic);
    }
    cout << "wall functions: calculate nodal normals at wall nodes." << endl;


    cout << "wall functions: set wall-function normals equal to nearest wall-node normals..." << endl;
    for (int inwf=0; inwf<(int) WFnodes.size(); inwf++) {
      WFnodes[inwf].modn = 0.;
      for (int id=0; id<Ndim; id++) {
        WFnodes[inwf].n[id]=No_norm[id][WFnodes[inwf].twin];
        WFnodes[inwf].modn += WFnodes[inwf].n[id]*WFnodes[inwf].n[id];
      }
      WFnodes[inwf].modn = sqrt(WFnodes[inwf].modn);
    }
    cout << "wall functions: set wall-function normals equal to nearest wall-node normals." << endl;

    cout << "wall functions initialization." << endl;
  }


  if (periodic) {
    cout << "find twins for periodic nodes..." << endl;
    for (int ig=1; ig<=Nbcgroup; ig++) {
      if (BCgroup[ig].type==IBPERI) {
        for (int inb=1; inb<=BCgroup[ig].nnode; inb++) {
          const int inu = Nobg[ig][inb].node;
          Nobg[ig][inb].twin=-1;
          for (int jg=1; jg<=Nbcgroup; jg++) {
            if (BCgroup[jg].type==IBPERE)
              for (int jnb=1; jnb<=BCgroup[jg].nnode; jnb++) {
                const int jnu = Nobg[jg][jnb].node;
                int flag = 0;
                for (int id=1; id<Ndim; id++) {
                  const int jd = (id+BCgroup[ig].option-1)%Ndim;
                  if (std::abs(M.vv[jd][inu]-M.vv[jd][jnu])<1.e-10)
                    flag++;
                }
                if (flag==Ndim-1) {
                  Nobg[ig][inb].twin=jnu;
                  break;
                }
              }
          }
          if (Nobg[ig][inb].twin==-1)
            nrerror("No twin found for periodic node");
        }
      }
      else if (BCgroup[ig].type==IBPERE) {
        for (int inb=1; inb<=BCgroup[ig].nnode; inb++) {
          const int inu = Nobg[ig][inb].node;
          Nobg[ig][inb].twin=-1;
          for (int jg=1; jg<=Nbcgroup; jg++) {
            if (BCgroup[jg].type==IBPERI)
              for (int jnb=1; jnb<=BCgroup[jg].nnode; jnb++) {
                const int jnu = Nobg[jg][jnb].node;
                int flag = 0;
                for (int id=1; id<Ndim; id++) {
                  const int jd = (id+BCgroup[ig].option-1)%Ndim;
                  if (std::abs(M.vv[jd][inu]-M.vv[jd][jnu])<1.e-10)
                    flag++;
                }
                if (flag==Ndim-1) {
                  Nobg[ig][inb].twin=jnu;
                  break;
                }
              }
          }
          if (Nobg[ig][inb].twin==-1)
            nrerror("No twin found for periodic node");
        }
      }
    }
    cout << "find twins for periodic nodes..." << endl;


    cout << "correcting periodic cells volumes..." << endl;
    for (int ig=1; ig<=Nbcgroup; ig++) {
      if (BCgroup[ig].type==IBPERI) {
        for (int inb=1; inb<=BCgroup[ig].nnode; inb++) {
          const int i = Nobg[ig][inb].node;
          const int j = Nobg[ig][inb].twin;
          No_vol[i] += No_vol[j];
          No_vol[j]  = No_vol[i];
        }
        break;
      }
    }
    cout << "correcting periodic cells volumes." << endl;


    cout << "set periodic gradient..." << endl;
    double pinlet = 0., xinlet = 0.,
           pexit  = 0., xexit  = 0.;
    for (int g=1; g<=Nbcgroup; g++) {
      if (BCgroup[g].type==IBPERI) {
        pinlet = BCgroup[g].invals[0];
        xinlet = M.vv[periodic_dirn][Nobg[g][1].node];
      }
      else if (BCgroup[g].type==IBPERE) {
        pexit = BCgroup[g].invals[0];
        xexit = M.vv[periodic_dirn][Nobg[g][1].node];
      }
    }
    periodic_pgrad = (pexit-pinlet)/(xexit-xinlet);
    cout << "set periodic gradient." << endl;
  }
}

