
/* read start file for input settings */

#include "common.h"


/*
 * read and setup from .c.xml testcase file, either from:
 * first <cdb><c/></cdb> node, or a specific <cdb><c label="label"/></cdb>
 * in: file[:label] with testcase
 */
void readstart(const std::string& ccase)
{
  using namespace std;


  // open file
  const string file  = ccase.substr(0,ccase.find(":"));
  const string label = ccase.find(":")==string::npos? "*" : ccase.substr(ccase.find(":")+1);
  cout << "readstart: reading " << file << ":<cdb><c label=\"" << label << "\"/></cdb>..." << endl;

  m::hal::xobject c = XMLNode::openFileHelper(file.c_str(),"cdb");
  c = label!="*"? c.getChildNodeWithAttribute("c","label",label.c_str())
                : c.getChildNode("c");
  if (c.isEmpty())
    nrerror("readstart: <cdb><c label=\"" + label + "\"/></cdb> not found");
  c.flatten();


  // <setup/>
  Ndim        = c.getChildNode("setup").getAttribute< int    >("dimensions");
  restart     = c.getChildNode("setup").getAttribute< string >("restart","false")=="true"? 1:0;
  file_input  = c.getChildNode("setup").getAttribute< string >("file_input");
  file_inlet  = c.getChildNode("setup").getAttribute< string >("file_inlet");
  file_output = c.getChildNode("setup").getAttribute< string >("file_output","c.plt");
  file_log = file_output + ".txt";


  // <terms/>
  rho         = c.getChildNode("terms").getAttribute< double >("density");
  nulam       = c.getChildNode("terms").getAttribute< double >("kviscosity");
  prandtl     = c.getChildNode("terms").getAttribute< double >("prandtl");
  grav[0]     = c.getChildNode("terms").getAttribute< double >("gx");
  grav[1]     = c.getChildNode("terms").getAttribute< double >("gy");
  grav[2]     = c.getChildNode("terms").getAttribute< double >("gz");
  diffusion   = c.getChildNode("terms").getAttribute< string >("diffusion",   "true")=="false"? 0:1;
  temperature = c.getChildNode("terms").getAttribute< string >("temperature", "false")=="true"? 1:0;
  buoyancy    = c.getChildNode("terms").getAttribute< string >("buoyancy",    "false")=="true"? 1:0;
  varden      = c.getChildNode("terms").getAttribute< string >("vardensity",  "false")=="true"? 1:0;
  rhofac      = c.getChildNode("terms").getAttribute< double >("vardensity_value");
  const string scalarscheme = c.getChildNode("terms").getAttribute< string >("scalarscheme","none");
  const string turbulence   = c.getChildNode("terms").getAttribute< string >("turbulence","none");
  scaconv = (scalarscheme=="FOU"?          ISSFOU :  // FOU finite-volume scheme
            (scalarscheme=="N"?            ISSNSC :  // N-scheme
            (scalarscheme=="Galerkin"?     ISSGAL :  // Galerkin
            (scalarscheme=="LDA"?          ISSLDA :  // LDA-scheme
            (scalarscheme=="Lax-Wendroff"? ISSLWS :  // Lax-Wendroff
            (scalarscheme=="PSI"?          ISSPSI :  // PSI-scheme
                                           ISNONE ))))));  // (no convection scheme)
  turmod = (turbulence=="KEHR"? ITKEHR :  // k-epsilon, high-Re with standard wall functions
           (turbulence=="KEHG"? ITKEHG :  // ... high-Re with generalized wall functions
           (turbulence=="KELB"? ITKELB :  // ... Lam-Bremhorst low-Re
           (turbulence=="KELS"? ITKELS :  // ... Launder-Sharma (Durbin version)
           (turbulence=="KENA"? ITKENA :  // ... Abe-Kondoh-Nagano low-Re
           (turbulence=="KE2L"? ITKE2L :  // ... two-layer
           (turbulence=="KWWF"? ITKWWF :  // k-omega, with wall functions
           (turbulence=="KWHR"? ITKWHR :  // ... Wilcox high-Re
           (turbulence=="KWLR"? ITKWLR :  // ... Wilcox low-Re
           (turbulence=="KWBS"? ITKWBS :  // ... Menter baseline model
           (turbulence=="KWSS"? ITKWSS :  // ... Menter shear-stress transport
           (turbulence=="KWPD"? ITKWPD :  // ... Peng-Davidson-Holmberg
           (turbulence=="LB"?   ITKELB :  // shortcuts for nice models
           (turbulence=="AKN"?  ITKENA :  // ...
           (turbulence=="WHR"?  ITKWHR :  // ...
           (turbulence=="WLR"?  ITKWLR :  // ...
           (turbulence=="BSL"?  ITKWBS :  // ...
           (turbulence=="SST"?  ITKWSS :  // ...
           (turbulence=="PDH"?  ITKWPD :  // ...
                                ITNONE )))))))))))))))))));  // (laminar)
  if (buoyancy && !temperature)
    nrerror("readstart: buoyancy terms only possible with temperature equation");


  // set number of cell/face nodes and equations, and wall distance properties
  Nvtcell  = 1 + Ndim;
  Nvtfce   = Ndim;
  dNvtcell = (double) Nvtcell;
  dNvtfce  = (double) Nvtfce;

  Ncoupled = 1 + Ndim;
  Neqns    = Ncoupled + (temperature? 1:0) + (turmod? 2:0);
  Nsys     = Neqns;  // temporary setting?
  iv_temp  = (!temperature? 0 : Ncoupled);
  iv_turb1 = (!turmod?      0 : Ncoupled + (temperature? 1:0));
  iv_turb2 = (!turmod?      0 : Ncoupled + (temperature? 1:0) + 1);

  wall_functions = turmod==ITKEHR || turmod==ITKEHG || turmod==ITKWWF? 1:0;
  walldist       = turmod && !wall_functions?                          1:0;


  // <mitremassembler/>
  mitremassembler.x  = c.getChildNode("mitremassembler");
  mitremassembler.ok = !(mitremassembler.x.isEmpty());
  if (mitremassembler.ok) {
    cout << "readstart: mitremassembler setup..." << endl;

    mitremassembler.m_mitrem = new MITReM(
      mitremassembler.x.getChildNode("MITReM").getAttribute< std::string >("file"),
      mitremassembler.x.getChildNode("MITReM").getAttribute< std::string >("label") );

    XMLNode ema = mitremassembler.x.getChildNode("ElementMatrixAssembler");
    mitremassembler.m_assembler = new ElementMatrixAssembler(
      (Ndim==2? "2D":"3D"),
      mitremassembler.m_mitrem,
      ema.getAttribute< std::string >("convectionScheme"),
      ema.getAttribute< std::string >("diffusionScheme"),
      ema.getAttribute< std::string >("migrationScheme"),
      ema.getAttribute< std::string >("magneticScheme"),
      ema.getAttribute< std::string >("homReactionScheme"),
      ema.getAttribute< std::string >("electrostaticsScheme"),
      ema.getAttribute< std::string >("timeScheme"),
      ema.getAttribute< std::string >("elecReactionScheme"),
      ema.getAttribute< std::string >("gasReactionScheme"),
      ema.getAttribute< std::string >("is_bubble","false")=="true",
      ema.getAttribute< std::string >("charge_flux","true")=="true",
      ema.getAttribute< std::string >("swap_first_and_last_equations","true")=="true" );

    XMLNode ls = mitremassembler.x.getChildNode("ls");
    mitremassembler.ls = m::Create< LS >(ls.getAttribute< string >("type","ls_aztec"));
    for (int a=0; a<ls.nAttribute(); ++a)
      (mitremassembler.ls)->xml.updateAttribute(ls.getAttribute(a).lpszValue,NULL,ls.getAttribute(a).lpszName);

    mitremassembler.iterinit = mitremassembler.x.getChildNode("iterinit").getAttribute< int >("value",0);
    mitremassembler.iv       = Neqns;
    mitremassembler.Nions    = (int) (mitremassembler.m_mitrem)->getNIons();
    mitremassembler.linrelx  = mitremassembler.x.getChildNode("linrelx").getAttribute< double >("value",1.);

    // override: force surface gas fraction, k and cSat
    XMLNode ss = mitremassembler.x.getChildNode("supersaturation");
    mitremassembler.surfacegasfraction_min = 0.01;
    mitremassembler.surfacegasfraction_max = 0.01;
    if (!(ss.isEmpty())) {
      double tmin = ss.getAttribute< double >("surfacegasfraction_min",0.01),
             tmax = ss.getAttribute< double >("surfacegasfraction_max",0.01);
      tmin = std::max(0.,std::min(1.,std::min(tmax,tmin)));
      tmax = std::max(0.,std::min(1.,std::max(tmax,tmin)));
      mitremassembler.surfacegasfraction_min = tmin;
      mitremassembler.surfacegasfraction_max = tmax;

      const double force_k    = ss.getAttribute< double >("k",1.),
                   force_cSat = ss.getAttribute< double >("cSat",0.);
      cout << " supersaturation (all gasreactions): force"
           << " k="    << force_k
           << " cSat=" << force_cSat << endl;
      for (unsigned r=0; r<(mitremassembler.m_mitrem)->getNGasReactions(); ++r) {
        (mitremassembler.m_mitrem)->setGasReactionKinParam(r,0,force_k);
        (mitremassembler.m_mitrem)->setGasReactionKinParam(r,1,force_cSat);
      }
    }

    // set bulk concentrations
    mitremassembler.forcebulk = mitremassembler.x.getAttribute< std::string >("forcebulk","false")=="true";
    mitremassembler.bulk.resize(mitremassembler.Nions);
    for (int i=0; i<mitremassembler.Nions; ++i)
      mitremassembler.bulk[i] = (mitremassembler.m_mitrem)->getIonInletConcentration(i);
    (mitremassembler.m_mitrem)->calcEquilibrium(mitremassembler.bulk);
    for (int i=0; i<mitremassembler.Nions; ++i)
      cout << " bulk(" << (mitremassembler.m_mitrem)->getIonLabel(i) << "): " << mitremassembler.bulk[i] << endl;

    cout << "readstart: mitremassembler setup." << endl;
  }
  else {
    cout << "readstart: mitremassembler not active." << endl;
  }
  Nmit  = mitremassembler.ok? (mitremassembler.Nions+1) : 0;


  // <solver_nonlinear/>
  Nitermax            = c.getChildNode("solver_nonlinear").getAttribute< int    >("maxiterations",1);
  Nitermin            = c.getChildNode("solver_nonlinear").getAttribute< int    >("miniterations",0);
  conv_thresh         = c.getChildNode("solver_nonlinear").getAttribute< double >("convergence_level");
  newton_eps          = c.getChildNode("solver_nonlinear").getAttribute< double >("newton_eps");
  newton_thresh       = c.getChildNode("solver_nonlinear").getAttribute< double >("newton_switch");
  turb_iterinit       = c.getChildNode("solver_nonlinear").getAttribute< int    >("turb_iterinit",5);
  scalar_coupling     = c.getChildNode("solver_nonlinear").getAttribute< string >("couple_t", "false")=="true"? 1:0;
  turbulence_coupling = c.getChildNode("solver_nonlinear").getAttribute< string >("couple_ke","false")=="true"? 1:0;
  const string convergence_variable = c.getChildNode("solver_nonlinear").getAttribute< string >("convergence_variable","p");
  const string method               = c.getChildNode("solver_nonlinear").getAttribute< string >("method","Newton");
  const string relax_puvw_type      = c.getChildNode("solver_nonlinear").getAttribute< string >("relax_puvw_type","linear");
  const string relax_t_type         = c.getChildNode("solver_nonlinear").getAttribute< string >("relax_t_type",   "linear");
  const string relax_ke_type        = c.getChildNode("solver_nonlinear").getAttribute< string >("relax_ke_type",  "linear");
  const double relax_puvw_value     = c.getChildNode("solver_nonlinear").getAttribute< double >("relax_puvw_value",1.);
  const double relax_t_value        = c.getChildNode("solver_nonlinear").getAttribute< double >("relax_t_value",   1.);
  const double relax_ke_value       = c.getChildNode("solver_nonlinear").getAttribute< double >("relax_ke_value",  1.);

  // correct real turbulent coefficients calculation
  turb_iterinit = restart? 0 : turb_iterinit;

  // set non-linear method
  Jacobian = (method=="Picard"?  0 :  // Picard
             (method=="Approx"?  1 :  // Newton, approximate
             (method=="Newton"?  2 :  // Newton (full)
                                -1 )));
  if (Jacobian<0)
    cout << "readstart: <solver_nonlinear method!=\"(Picard|Approx|Newton)\" />, coupled system calculation will be skipped!" << endl;

  // set linear or global/local time-step relaxation
  dtrelax        = (relax_puvw_type=="dt-global"? 1 : (relax_puvw_type=="dt-local"?  2 : 0 ));
  dtrelax_scalar = (relax_t_type   =="dt-global"? 1 : (relax_t_type   =="dt-local"?  2 : 0 ));
  dtrelax_turb   = (relax_ke_type  =="dt-global"? 1 : (relax_ke_type  =="dt-local"?  2 : 0 ));
  CFL            = !dtrelax?        0. : relax_puvw_value;
  CFL_scalar     = !dtrelax_scalar? 0. : relax_t_value;
  CFL_turb       = !dtrelax_turb?   0. : relax_ke_value;
  linrlx         =  dtrelax?        1. : relax_puvw_value;
  linrlx_scalar  =  dtrelax_scalar? 1. : relax_t_value;
  linrlx_turb    =  dtrelax_turb?   1. : relax_ke_value;

  // correct number of equations based on couplings
  if (temperature && scalar_coupling)      Ncoupled += 1;
  if (turmod && (turbulence_coupling==2))  Ncoupled += 2;


  // <system_(coupled|scalar|turb)/>
  if (temperature && !scalar_coupling) {
    XMLNode ls = c.getChildNode("system_scalar");
    if (ls.isEmpty())
      nrerror("readstart: no system_scalar present, but it is required");
    ls_scalar = m::Create< LS >(ls.getAttribute< string >("type","ls_aztec"));
    for (int a=0; a<ls.nAttribute(); ++a)
      ls_scalar->xml.updateAttribute(ls.getAttribute(a).lpszValue,NULL,ls.getAttribute(a).lpszName);
  }
  if (turmod && turbulence_coupling==0) {
    XMLNode ls = c.getChildNode("system_turb1");
    if (ls.isEmpty())
      nrerror("readstart: no system_turb1 present, but it is required");
    ls_turb1 = m::Create< LS >(ls.getAttribute< string >("type","ls_aztec"));
    for (int a=0; a<ls.nAttribute(); ++a)
      ls_turb1->xml.updateAttribute(ls.getAttribute(a).lpszValue,NULL,ls.getAttribute(a).lpszName);
    ls = c.getChildNode("system_turb2");
    if (ls.isEmpty())
      nrerror("readstart: no system_turb2 present, but it is required");
    ls_turb2 = m::Create< LS >(ls.getAttribute< string >("type","ls_aztec"));
    for (int a=0; a<ls.nAttribute(); ++a)
      ls_turb2->xml.updateAttribute(ls.getAttribute(a).lpszValue,NULL,ls.getAttribute(a).lpszName);
  }
  if (turmod && turbulence_coupling==1) {
    XMLNode ls = c.getChildNode("system_turb");
    if (ls.isEmpty())
      nrerror("readstart: no system_turb present, but it is required");
    ls_turb = m::Create< LS >(ls.getAttribute< string >("type","ls_aztec"));
    for (int a=0; a<ls.nAttribute(); ++a)
      ls_turb->xml.updateAttribute(ls.getAttribute(a).lpszValue,NULL,ls.getAttribute(a).lpszName);
  }
  if (true) {
    XMLNode ls = c.getChildNode("system_coupled");
    if (ls.isEmpty())
      nrerror("readstart: no system_coupled present, but it is required");
    ls_coupled = m::Create< LS >(ls.getAttribute< string >("type","ls_aztec"));
    for (int a=0; a<ls.nAttribute(); ++a)
      ls_coupled->xml.updateAttribute(ls.getAttribute(a).lpszValue,NULL,ls.getAttribute(a).lpszName);
  }


  // <vars/>
  if (c.getChildNode("vars").nChildNode("v")!=Neqns)
    nrerror("readstart: number of variables not correct with the settings");
  m_vars_label.assign(Neqns+Nmit,"");
  m_vars_init .assign(Neqns+Nmit,0.);
  for (int i=0; i<Neqns; ++i) {
    XMLNode v = c.getChildNode("vars").getChildNode("v",i);
    m_vars_label[i] = v.getAttribute< string >("label","?");
    m_vars_init [i] = v.getAttribute< double >("init",0.);
  }
  if (Nmit) {
    for (int i=0; i<mitremassembler.Nions; ++i) {
      m_vars_label[Neqns+i] = (mitremassembler.m_mitrem)->getIonLabel(i);
      m_vars_init [Neqns+i] = (mitremassembler.m_mitrem)->getIonInletConcentration(i);
    }
    m_vars_label[Neqns+mitremassembler.Nions] = "UU";
    m_vars_init [Neqns+mitremassembler.Nions] = 0.;
  }
  logL1 = logL2 = logLi = std::vector< double >(Neqns+Nmit,0.);
  resL1 = resL2 = resLi = std::vector< double >(Neqns+Nmit,0.);


  // set convergence variable
  iverr = 0;
  for (int i=0; i<(int) m_vars_label.size(); ++i)
    if (convergence_variable==m_vars_label[i]) {
      iverr = i;
      break;
    }
  cout << "readstart: convergence variable: \"" << m_vars_label[iverr] << '"' << endl;


  // <bcs/>
  // resize vector first with default options (first comes a dummy)
  Nbcgroup = c.getChildNode("bcs").nChildNode("bc");
  BCgroup.resize(1+Nbcgroup);
  BCgroup.front().type = IBNONE;

  // set boundaries
  int ib = 0;
  for (vector< bcgroup_struct >::iterator b=BCgroup.begin()+1; b!=BCgroup.end(); ++b, ++ib) {

    // set type and its options
    b->zone             = c.getChildNode("bcs").getChildNode("bc",ib).getAttribute< string >("zone");
    const string type   = c.getChildNode("bcs").getChildNode("bc",ib).getAttribute< string >("type","none");
    const string option = c.getChildNode("bcs").getChildNode("bc",ib).getAttribute< string >("option","none");
    b->type = (type=="fixv"? IBFIXV :
              (type=="fixp"? IBFIXP :
              (type=="wall"? IBWALL :
              (type=="symm"? IBSYMM :
              (type=="peri"? IBPERI :
              (type=="pere"? IBPERE :
                             IBNONE ))))));
    switch (b->type) {
      case (IBFIXV): {
        if      (option=="u") { b->option=0;  b->invals.resize(1 + (temperature? 1:0) +(turmod? 2:0)); }  // uniform
        else if (option=="p") { b->option=1;  b->invals.resize(3 + (temperature? 1:0) +(turmod? 2:0)); }  // parabolic
        else if (option=="b") { b->option=2;  b->invals.resize(3 + (temperature? 1:0) +(turmod? 2:0)); }  // boundary layer
        else if (option=="r") { b->option=3;  b->invals.resize(Ndim>2? 5:2);                           }  // read from file
        else if (option=="c") { b->option=4;  b->invals.resize(5 + (temperature? 1:0) +(turmod? 2:0)); }  // parabolic circular-duct
        else
          nrerror("readstart: fixv unknown option");
        if (option=="c" && Ndim!=3)
          nrerror("readstart: fixv parabolic circular-duct only for 3D");
      } break;
      case (IBWALL): {
        b->type   = (option=="q") && temperature? IBWALQ:IBWALL;  // heat flux wall?
        b->option =  option=="m"? 1:0;                            // moving wall?
        b->invals.resize((b->option? Ndim:0) + (temperature? 1:0));
      } break;
      case (IBFIXP): {
        b->invals.resize(1.);
        if      (option=="x")  b->option=1;
        else if (option=="y")  b->option=2;
        else if (option=="z")  b->option=3;
        else                   b->option=0;
      } break;
      case (IBSYMM): {
        b->invals.clear();
        if      (option=="x")  b->option=1;
        else if (option=="y")  b->option=2;
        else if (option=="z")  b->option=3;
        else
          nrerror("readstart: symm unknown option");
      } break;
      case (IBPERI): case (IBPERE): {
        b->invals.resize(1.);
        if      (option=="x")  b->option=1;
        else if (option=="y")  b->option=2;
        else if (option=="z")  b->option=3;
        else
          nrerror("readstart: peri/pere unknown option");
      } break;
      default: break;
    }

    // set values
    if (b->invals.size()) {
      istringstream values( c.getChildNode("bcs").getChildNode("bc",ib).getAttribute("values") );
      for (unsigned i=0; i<b->invals.size(); ++i) {
        values >> b->invals[i];
        if (!values)
          nrerror("readstart: boundary has incorrect number of values");
      }
    }
  }

  // add an extra entry for wall functions and set periodicity
  if (wall_functions) {
    BCgroup.push_back(BCgroup.front());
    BCgroup.back().type = IBWLFN;
  }
  periodic = 0;
  for (vector< bcgroup_struct >::iterator b=BCgroup.begin(); b!=BCgroup.end(); ++b)
    if (b->type==IBPERI) {
      periodic      = 1;
      periodic_dirn = b->option - 1;
      break;
    }


  // initialize parameters and turbulence model constants
  To    = 0.;
  Qin   = 0.;
  Qout  = 0.;
  iter  = 0;
  epsilon = 1.e-20;
  if (turmod)
    turb_init((ITid) turmod);


  cout << "readstart: reading." << endl;
}

