
#include <numeric>
#include <memory>
#include <algorithm>
#include "boost/progress.hpp"
#include "common.h"


// forward declarations
enum CALL { ASSEMBLY, CURRENT, GAS, SIZE, BSIZE };
void assembleMITReM(const std::vector< m::melem >& e2n);
void assembleDirichlet(const std::vector< m::melem >& e2n,
  const std::vector< unsigned >& vlabels,
  const std::vector< double   >& vvalues);
double assembleElectrode(
  const std::vector< m::melem >& e2n,
  const std::vector< unsigned >& greactions,
  const std::vector< unsigned >& ereactions,
  const double& v,
  const CALL& action=ASSEMBLY);

// forward declarations (utility functions)
std::vector< unsigned > getZones(const std::string& l);
std::vector< unsigned > getVLabels(const std::string& l);
std::vector< double   > getVValues(const std::string& l);
std::vector< unsigned > getEReactions(const std::string& l);
std::vector< unsigned > getGReactions(const std::string& l);
std::vector< std::string > getVStrings(const std::string& l);


using namespace std;


void Update_mitrem()
{
  // - zero jacobian matrix, solution and residual vectors
  // - set main connectivity in the right place...
  // - create vector of indices for (all) bulk concentrations
  // - create vector of indices for (all) electrode reactions
  mitremassembler_struct& m = mitremassembler;
  (m.ls)->reset();
  e2n.swap(M.vz[0].e2n);
  vector< unsigned > ibulk(m.Nions,0);
  for (int i=1; i<m.Nions; ++i)
    ibulk[i] = ibulk[i-1]+1;


  {
    boost::progress_timer time(cout);
    cout << "Update_mitrem: system assembly..." << endl;


    cout << "Update_mitrem: assemble MITReM equations..." << endl;
    assembleMITReM(M.vz[0].e2n);
    cout << "Update_mitrem: assemble MITReM equations." << endl;


    for (int i=0; i<m.x.getChildNode("bcs").nChildNode("bc"); ++i) {
      const XMLNode b = m.x.getChildNode("bcs").getChildNode("bc",i);
      const string  t = b.getAttribute("type");
      const string  l = b.getAttribute("label");
      const vector< unsigned > z = getZones(b.getAttribute< string >("zone"));
      for (vector< unsigned >::const_iterator j=z.begin(); j!=z.end(); ++j) {
        cout << "Update_mitrem: assemble bc \"" << l << "\" (" << t << ") zone \"" << M.vz[*j].n << "\"..." << endl;

        if      (t=="bulk")      assembleDirichlet( M.vz[*j].e2n,ibulk,m.bulk );
        else if (t=="velectrode") {
          // vector of variable indices, and values
          vector< unsigned > ibulk_and_u(ibulk);
          vector< double   > bulk_and_u(m.bulk);
          ibulk_and_u.push_back(ibulk.back()+1);
          bulk_and_u.push_back(b.getAttribute< double >("metalpotential",0.));
          assembleDirichlet( M.vz[*j].e2n,ibulk_and_u,bulk_and_u );
        }
        else if (t=="dirichlet") assembleDirichlet( M.vz[*j].e2n,
          getVLabels(b.getAttribute< string >("vlabels")),
          getVValues(b.getAttribute< string >("vvalues")) );
        else if (t=="electrode") assembleElectrode( M.vz[*j].e2n,
          getGReactions(b.getAttribute< string >("gasreaction")),
          getEReactions(b.getAttribute< string >("elecreaction")),
          b.getAttribute< double >("metalpotential") );

        cout << "Update_mitrem: assemble bc \"" << l << "\" (" << t << ") zone \"" << M.vz[*j].n << "\"." << endl;
      }
    }


    cout << "Update_mitrem: system assembly: timer: ";
  }

  // update residual L2 error norm
  double residu = 0.;
  for (unsigned c=0; c<(m.ls)->Nb*(m.ls)->Ne; ++c)
    residu += (m.ls)->B(c)*(m.ls)->B(c);
  residu = sqrt(residu);
  cout << "Update_mitrem: residual = " << residu << endl;


  // set back main connectivity
  e2n.swap(M.vz[0].e2n);


  cout << "Update_mitrem: solve linear system..." << endl;
  {
    boost::progress_timer time(cout);
    (m.ls)->solve();
    cout << "Update_mitrem: timer: ";
  }
  cout << "Update_mitrem: solve linear system." << endl;


  cout << "Update_mitrem: update solution..." << endl;
  // update solution, count neg. concentrations and solution L2 error norm
  vector< unsigned > negc(m.Nions,0);
  for (int i=0; i<m.Nions+1; ++i) {
    double l2 = 0.;
    for (int n=0; n<Nnode; ++n) {
      const double dx = (m.ls)->X(n,i);
      No_W[m.iv+i][n] += dx*m.linrelx;
      l2              += dx*dx;
      if (i<m.Nions && No_W[m.iv+i][n]<0.)  ++negc[i];
    }
    logL2[m.iv+i] = max( log10(sqrt(l2/(double) Nnode)), -42. );
  }
  for (int g=1; periodic && g<=Nbcgroup; ++g) {
    if (BCgroup[g].type==IBPERE)
      for (int inb=1; inb<=BCgroup[g].nnode; ++inb)
        for (int iv=m.iv; iv<m.Nions+1; ++iv)
          No_W[iv][Nobg[g][inb].node] = No_W[iv][Nobg[g][inb].twin];
  }
  const unsigned negc_sum = accumulate(negc.begin(),negc.end(),0);
  if (negc_sum) {
    cout << "Update_mitrem: neg. concentrations sum: " << negc_sum << endl;
    for (int i=0; i<m.Nions; ++i)
      if (negc[i])
        cout << "Update_mitrem: neg. concentrations (\"" << (m.m_mitrem)->getIonLabel((unsigned) i) << "\"): " << negc[i] << endl;
  }
  cout << "Update_mitrem: update solution." << endl;


  // calculate current, current densities and gas production rate
  for (unsigned i=0; i<(m.m_mitrem)->getNElecReactions(); ++i)
    m.Mj.vv[Ndim+i].assign(Nnode,0.);
  for (unsigned i=0; i<(m.m_mitrem)->getNElecReactions()+(m.m_mitrem)->getNGasReactions(); ++i)
    m.Mv.vv[Ndim+i].assign(Nnode,0.);
  {
    // create per-zone-indexed vectors to store metal potentials and
    // electrode/gas reaction indices (removing duplicate entries)
    vector< double >              vz_mp(M.z(),0.);
    vector< vector < unsigned > > vz_er(M.z()),
                                  vz_gr(M.z());
    for (int i=0; i<m.x.getChildNode("bcs").nChildNode("bc"); ++i) {
      const XMLNode b = m.x.getChildNode("bcs").getChildNode("bc",i);
      const string  t = b.getAttribute("type");
      if (t=="electrode") {
        const vector< unsigned > z = getZones(b.getAttribute< string >("zone"));
        for (vector< unsigned >::const_iterator j=z.begin(); j!=z.end(); ++j) {
          vector< unsigned > er(getEReactions(b.getAttribute< string >("elecreaction"))),
                             gr(getGReactions(b.getAttribute< string >("gasreaction")));
          vz_mp[*j] = b.getAttribute< double >("metalpotential");
          vz_er[*j].insert(vz_er[*j].end(),er.begin(),er.end());
          vz_gr[*j].insert(vz_gr[*j].end(),gr.begin(),gr.end());
        }
      }
    }
    for (unsigned i=1; i<M.z(); ++i) {
      sort(vz_er[i].begin(),vz_er[i].end());
      sort(vz_gr[i].begin(),vz_gr[i].end());
      vz_er[i].erase(unique(vz_er[i].begin(),vz_er[i].end()),vz_er[i].end());
      vz_gr[i].erase(unique(vz_gr[i].begin(),vz_gr[i].end()),vz_gr[i].end());
    }

    cout << "Update_mitrem: calculating electrochemical properties..." << endl;
    double Ibalance = 0;
    for (unsigned i=1; i<M.z(); ++i) {
      if (!vz_gr[i].size() && !vz_er[i].size())
        continue;
      const double I = assembleElectrode(M.vz[i].e2n,vz_gr[i],vz_er[i],vz_mp[i],CURRENT),
                   J = I/max(assembleElectrode(M.vz[i].e2n,vz_gr[i],vz_er[i],vz_mp[i],BSIZE),1.e-20),
                   G = assembleElectrode(M.vz[i].e2n,vz_gr[i],vz_er[i],vz_mp[i],GAS);
      Ibalance += I;
      cout << "Update_mitrem: zone: \"" << M.vz[i].n << '"'
           << "  I[A]="      << I
           << "  J[A.m-2]="  << J
           << "  G[m3.s-1]=" << G << endl;
    }
    cout << "Update_mitrem: current balance [A]=" << Ibalance << endl;
    cout << "Update_mitrem: calculating electrochemical properties." << endl;
  }


  {
    string outfile_j = file_output + ".j" + extension(file_output),
           outfile_v = file_output + ".v" + extension(file_output);


    cout << "Update_mitrem: writing current densities to \"" << outfile_j << "\"..." << endl;
    m.Mj.vz.swap(M.vz);          // (hack)
    m.Mj.vz[0].e2n.swap(e2n);    // ...
    for (int d=0; d<Ndim; ++d)   // ...
      m.Mj.vv[d].swap(M.vv[d]);  // ...

    auto_ptr< m::mfoutput > p(m::mfactory< m::mfoutput  >::instance()->Create(extension(outfile_j)));
    {
      char* argv[] = { (char*) "", const_cast< char* >(outfile_j.c_str()) };
      GetPot o2(2,argv);
      p->write(o2,m.Mj);
    }

    for (int d=0; d<Ndim; ++d)   // (hack back)
      m.Mj.vv[d].swap(M.vv[d]);  // ...
    m.Mj.vz[0].e2n.swap(e2n);    // ...
    m.Mj.vz.swap(M.vz);          // ...
    cout << "Update_mitrem: writing current densities." << endl;


    cout << "Update_mitrem: writing electrode/gas reaction rates to \"" << outfile_v << "\"..." << endl;
    m.Mv.vz.swap(M.vz);          // (hack again)
    m.Mv.vz[0].e2n.swap(e2n);    // ...
    for (int d=0; d<Ndim; ++d)   // ...
      m.Mv.vv[d].swap(M.vv[d]);  // ...

    {
      char* argv[] = { (char*) "", const_cast< char* >(outfile_v.c_str()) };
      GetPot o2(2,argv);
      p->write(o2,m.Mv);
    }

    for (int d=0; d<Ndim; ++d)   // (hack back)
      m.Mv.vv[d].swap(M.vv[d]);  // ...
    m.Mv.vz[0].e2n.swap(e2n);    // ...
    m.Mv.vz.swap(M.vz);          // ...
    cout << "Update_mitrem: writing electrode/gas reaction rates." << endl;
  }
}


void assembleMITReM(const vector< m::melem >& e2n_)
{
  // some shortcuts
  mitremassembler_struct& m = mitremassembler;
  const int Nenod = (int) e2n_[0].n.size();
  const int Nions = m.Nions;
  const int Nsys_ = m.Nions+1;
  const vector< double > temperatures (Nenod,(m.m_mitrem)->getSolutionTemperature());
  const vector< double > densities    (Nenod,(m.m_mitrem)->getSolutionDensity());
  const vector< double > voidfractions(Nenod,0.);

  // allocate auxiliary variables
  double **coordinates    = dmatrix(0,Nenod-1,0,Ndim-1);
  double **velocities     = dmatrix(0,Nenod-1,0,Ndim-1);
  double **concentrations = dmatrix(0,Nenod-1,0,Nions-1);
  double **bvectors       = dmatrix(0,Nenod-1,0,3-1);
  vector< double > potentials(Nenod,0.);

  // element matrix and residual vector contribution
  double **emat = dmatrix(0,Nenod*Nsys_-1,0,Nenod*Nsys_-1);
  vector< double > eres(Nenod*Nsys_,0.);

  // cycle elements
  boost::progress_display pbar((unsigned) e2n_.size());
  for (vector< m::melem >::const_iterator e=e2n_.begin(); e!=e2n_.end(); ++e, ++pbar) {

    // set auxiliary variables
    for (int n=0; n<Nenod; ++n) {
      const unsigned N = (e->n)[n];
      for (int d=0; d<Ndim;  ++d) coordinates   [n][d] = M.vv[d][N];
      for (int d=0; d<Ndim;  ++d) velocities    [n][d] = No_W[1+d]   [N];
      for (int i=0; i<Nions; ++i) concentrations[n][i] = No_W[m.iv+i][N];
      for (int d=0; d<3;     ++d) bvectors      [n][d] = 0.;
      potentials[n] = No_W[m.iv+Nions][N];
    }

    // assemble element: Ax = b
    // (get element matrix contribution)
    DoubleMatrix emat_ = (m.m_assembler)->calcElementMat(
      coordinates, velocities, concentrations, &potentials[0],
      &temperatures[0], &densities[0], &voidfractions[0], bvectors );

    // copy into local matrix
    for (int ir=0; ir<Nenod*Nsys_; ++ir)
      for (int ic=0; ic<Nenod*Nsys_; ++ic)
        emat[ir][ic] = emat_[ir][ic];

    // element residual vector contribution (-Ax part of r = b-Ax)
    for (int ir=0; ir<Nenod*Nsys_; ++ir) {
      double sum = 0.;
      for (int inc=0; inc<Nenod; ++inc) {
        for (int i=0; i<Nions; ++i)
          sum += emat_[ir][inc*Nsys_+i]*concentrations[inc][i];
        sum += emat_[ir][inc*Nsys_+Nions]*potentials[inc];
      }
      eres[ir] = sum;
    }

    // assemble element: K = AMat + AJac - dB/dX (K*dx = r)
    // (get element jacobian matrix contribution)
    emat_ = (m.m_assembler)->calcElementJac(
      coordinates, velocities, concentrations, &potentials[0],
      &temperatures[0], &densities[0], &voidfractions[0], bvectors );

    // add to element contribution
    for (int ir=0; ir<Nenod*Nsys_; ++ir)
      for (int ic=0; ic<Nenod*Nsys_; ++ic)
        emat[ir][ic] += emat_[ir][ic];

    // add to system matrix and residual vector
    for (int iR=0; iR<Nenod; ++iR)
      for (int ir=0; ir<Nsys_; ++ir) {
        for (int iC=0; iC<Nenod; ++iC)
          for (int ic=0; ic<Nsys_; ++ic)
            (m.ls)->A((e->n)[iR],(e->n)[iC],ir,ic) += emat[iR*Nsys_+ir][iC*Nsys_+ic];
        (m.ls)->B((e->n)[iR],ir) -= eres[iR*Nsys_+ir];
      }
  }

  // deallocate auxiliary variables
  free_dmatrix(emat          ,0,Nenod*Nsys_-1,0,Nenod*Nsys_-1);
  free_dmatrix(bvectors      ,0,Nenod-1,0,3-1);
  free_dmatrix(concentrations,0,Nenod-1,0,Nions-1);
  free_dmatrix(velocities    ,0,Nenod-1,0,Ndim-1);
  free_dmatrix(coordinates   ,0,Nenod-1,0,Ndim-1);
}


void assembleDirichlet(
  const vector<m::melem> &e2n_,
  const vector< unsigned > &vlabels,
  const vector< double   > &vvalues)
{
  mitremassembler_struct& m = mitremassembler;
  for (vector< m::melem >::const_iterator e=e2n_.begin(); e!=e2n_.end(); ++e)
    for (vector< unsigned >::const_iterator n=e->n.begin(); n!=e->n.end(); ++n) {
      // variable index and value iterators
      vector< unsigned >::const_iterator i = vlabels.begin();
      vector< double   >::const_iterator v = vvalues. begin();
      for (; i!=vlabels.end() && v!=vvalues.end(); ++i, ++v) {
        (m.ls)->zerorow(*n,*i);
        (m.ls)->A(*n,*n,*i,*i) = 1.;
        (m.ls)->B(*n,*i)       = *v - No_W[m.iv+*i][*n];
      }
    }
}


double assembleElectrode(
  const vector< m::melem >& e2n_,
  const vector< unsigned >& greactions,
  const vector< unsigned >& ereactions,
  const double& v,
  const CALL& action)
{
  double ret = 0.;

  // some shortcuts
  mitremassembler_struct& m = mitremassembler;
  const int Nenod = (int) e2n_[0].n.size();
  const int Nions = m.Nions;
  const int Nsys_ = m.Nions+1;
  const unsigned* p_ereactions = (ereactions.size()? &ereactions[0]:NULL);
  const unsigned* p_greactions = (greactions.size()? &greactions[0]:NULL);

  // allocate auxiliary variables
  const vector< double > temperatures(Nenod,(m.m_mitrem)->getSolutionTemperature());
  const vector< double > densities   (Nenod,(m.m_mitrem)->getSolutionDensity());
  double **coordinates    = dmatrix(0,Nenod-1,0,Ndim-1);
  double **concentrations = dmatrix(0,Nenod-1,0,Nions-1);
  vector< double > potentials   (Nenod);
  vector< double > sgasfractions(Nenod,0.);

  // cycle elements
  for (vector< m::melem >::const_iterator e=e2n_.begin(); e!=e2n_.end(); ++e) {

    // set auxiliary variables
    for (int n=0; n<Nenod; ++n) {
      const unsigned N = (e->n)[n];
      for (int d=0; d<Ndim;  ++d) coordinates   [n][d] = M.vv[d][N];
      for (int i=0; i<Nions; ++i) concentrations[n][i] = No_W[m.iv+i][N];
      potentials   [n] = No_W[m.iv+Nions][N];
      if (greactions.size())
        sgasfractions[n] = min(m.surfacegasfraction_max,
                           max(m.surfacegasfraction_min,sgasfractions[n]));
    }

    if (action==ASSEMBLY) {
      // assemble: get face matrix and residual vector contributions
      DoubleMatrix fmat = (m.m_assembler)->calcBoundaryElementJac(
        coordinates, concentrations, &potentials[0],
        &temperatures[0], &densities[0], &sgasfractions[0],
        ereactions.size(), p_ereactions, v,
        greactions.size(), p_greactions );
      DoubleVector fres = (m.m_assembler)->calcBoundaryElementVec(
        coordinates, concentrations, &potentials[0],
        &temperatures[0], &densities[0], &sgasfractions[0],
        ereactions.size(), p_ereactions, v,
        greactions.size(), p_greactions );

      // add to system matrix and residual vector (b part of r = b-Ax)
      for (int iR=0; iR<Nenod; ++iR)
        for (int ir=0; ir<Nsys_; ++ir) {
          for (int iC=0; iC<Nenod; ++iC)
            for (int ic=0; ic<Nsys_; ++ic)
              (m.ls)->A((e->n)[iR],(e->n)[iC],ir,ic) += fmat[iR*Nsys_+ir][iC*Nsys_+ic];
          (m.ls)->B((e->n)[iR],ir) += fres[iR*Nsys_+ir];
        }
    }
    else if (action==CURRENT) {
      // calculate total current (by accumulation)
      ret += (m.m_assembler)->calcCurrent(
        coordinates, concentrations, &potentials[0],
        &temperatures[0], &densities[0], &sgasfractions[0],
        ereactions.size(), p_ereactions, v );

      // calculate current density and electrode/gas reaction rates
      DoubleListList j = (m.m_assembler)->calcElecReactionCurrentDensities(
        coordinates, concentrations, &potentials[0],
        &temperatures[0], &densities[0], &sgasfractions[0],
        ereactions.size(), p_ereactions, v );
      for (int i=0; i<Nenod; ++i) {
        for (unsigned r=0; r<ereactions.size(); ++r)
          m.Mj.vv[Ndim+ereactions[r]][(e->n)[i]] = j[i][r];
        (m.m_mitrem)->init(concentrations[i],potentials[i],temperatures[i],densities[i]);
        for (unsigned r=0; r<ereactions.size(); ++r)
          m.Mv.vv[Ndim+ereactions[r]][(e->n)[i]] =
            (m.m_mitrem)->calcElecReactionRate(r,v);
        for (unsigned r=0; r<greactions.size(); ++r)
          m.Mv.vv[Ndim+(m.m_mitrem)->getNElecReactions()+greactions[r]][(e->n)[i]] =
            (m.m_mitrem)->calcGasReactionRate(r);
      }
    }
    else if (action==GAS) {
      // calculate total gas production rate (by accumulation)
      ret += (m.m_assembler)->calcGasGeneration(
        coordinates, concentrations, &potentials[0],
        &temperatures[0], &densities[0], &sgasfractions[0],
        greactions.size(), p_greactions);
    }
    else if (action==SIZE) {
      // calculate total volume/area (thus size) by accumulation
      ret += (m.m_assembler)->calcESize(coordinates);
    }
    else if (action==BSIZE) {
      // calculate total area/length (thus size) by accumulation
      ret += (m.m_assembler)->calcBESize(coordinates);
    }
  }

  // deallocate auxiliary variables
  free_dmatrix(concentrations,0,Nenod-1,0,Nions-1);
  free_dmatrix(coordinates   ,0,Nenod-1,0,Ndim-1);

  return ret;
}


vector< unsigned > getZones(const string& l)
{
  vector< unsigned > r;
  vector< string > list = getVStrings(l);
  for (unsigned i=0; i<M.z(); ++i)
    for (unsigned j=0; j<(unsigned) list.size(); ++j)
      if (M.vz[i].n==list[j]) {
        r.push_back(i);
        break;
      }
  return r;
}


vector< unsigned > getVLabels(const string& l)
{
  vector< unsigned > r;
  vector< string > list = getVStrings(l);
  for (vector< string >::const_iterator s=list.begin(); s!=list.end(); ++s)
    if (*s=="UU")  r.push_back((unsigned) mitremassembler.Nions);
    else {
      for (unsigned i=0; i<(mitremassembler.m_mitrem)->getNIons(); ++i)
        if ((mitremassembler.m_mitrem)->getIonLabel(i)==*s) {
          r.push_back(i);
          break;
        }
    }
  return r;
}


vector< double > getVValues(const string& l)
{
  vector< double > r;
  vector< string > list = getVStrings(l);
  for (vector< string >::const_iterator s=list.begin(); s!=list.end(); ++s) {
    istringstream iss(*s);
    r.push_back(0.);
    iss >> r.back();
  }
  return r;
}


vector< unsigned > getEReactions(const string& l)
{
  vector< unsigned > r;
  vector< string > list = getVStrings(l);
  for (unsigned i=0; i<(mitremassembler.m_mitrem)->getNElecReactions(); ++i)
    for (unsigned j=0; j<(unsigned) list.size(); ++j)
      if ((mitremassembler.m_mitrem)->getElecReactionLabel(i)==list[j]) {
        r.push_back(i);
        break;
      }
  return r;
}


vector< unsigned > getGReactions(const string& l)
{
  vector< unsigned > r;
  vector< string > list = getVStrings(l);
  for (unsigned i=0; i<(mitremassembler.m_mitrem)->getNGasReactions(); ++i)
    for (unsigned j=0; j<(unsigned) list.size(); ++j)
      if ((mitremassembler.m_mitrem)->getGasReactionLabel(i)==list[j]) {
        r.push_back(i);
        break;
      }
  return r;
}


vector< string > getVStrings(const string& l)
{
  // split string by ':'
  string::size_type p1 = 0,
                    p2 = 0;
  vector< string > r;
  while (p2!=string::npos) {
    p2 = l.find(":",p1);
    istringstream ss( l.substr(p1,(p2==string::npos? p2:p2-p1)) );
    r.push_back(string());
    ss >> r.back();
    p1 = p2+1;
  }
  return r;
}
