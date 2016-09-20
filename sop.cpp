#include <iostream>
#include <fstream>
#include <cstring>
#include <math.h>
#include <cstdlib>
#include <ctime>
#include <cstdio>
#include <unistd.h>
#include "sop.h"
#include "random_generator.h"

// #define UNCORRECTED
const int debug = 0;

int ncmd;
char cmd[mcmd+1][mwdsize];
char opt[mopt_tot+1][mwdsize];
int opt_ptr[mcmd+1];
char pathname[MAXPATHLEN];

// bonded info

double k_bnd; // bond spring constant
int nbnd; // number of bonds
int* ibead_bnd;
int* jbead_bnd;
double* pdb_dist; // pdb bond distances
int bnds_allocated = 0;
double R0;
double R0sq;
double e_bnd_coeff;

// angular info

double k_ang;
int nang;
int* ibead_ang;
int* jbead_ang;
int* kbead_ang;
double* pdb_ang;
int angs_allocated = 0;
double e_ang_coeff;
double e_ang_ss_coeff;
double f_ang_ss_coeff;

// rna-rna vdw

int ncon_att; // number of native contacts
int ncon_rep; // repulisve non-native contact

// neighbor list
int nnl_att;
int nnl_rep;

// pair list
int nil_att;
int nil_rep;


double coeff_att[3][3] = { {0.0, 0.0, 0.0},
			   {0.0, 0.7, 0.8},
			   {0.0, 0.8, 1.0} };

double coeff_rep[3][3] = { {0.0, 0.0, 0.0},
			   {0.0, 1.0, 1.0},
			   {0.0, 1.0, 1.0} };

double force_coeff_att[3][3] = { {0.0,       0.0,       0.0},
				 {0.0, -12.0*1.0, -12.0*0.8},
				 {0.0, -12.0*0.8, -12.0*0.7} };

double force_coeff_rep[3][3] = { {0.0,       0.0,       0.0},
				 {0.0,  -6.0*1.0,  -6.0*1.0},
				 {0.0,  -6.0*1.0,  -6.0*1.0} };

double sigma_rep[3][3] = { {0.0, 0.0, 0.0},
			   {0.0, 3.8, 5.4},
			   {0.0, 5.4, 7.0} };

double sigma_rep2[3][3];
double sigma_rep6[3][3];
double sigma_rep12[3][3];

//double sigma_rep;
//double sigma_rep2;
double sigma_ss; // for angular soft-sphere repulsion
double sigma_ss6; // for angular soft-sphere repulsion
double epsilon_ss; // for angular soft-sphere repulsion
//double force_coeff_rep;
double rcut_nat[3][3] = { { 0.0,  0.0,  0.0},
                          { 0.0,  8.0, 11.0},
                          { 0.0, 11.0, 14.0} };
int* ibead_lj_nat;
int* jbead_lj_nat;
int* itype_lj_nat;
int* jtype_lj_nat;

double* lj_nat_pdb_dist;
double* lj_nat_pdb_dist2;
double* lj_nat_pdb_dist6;
double* lj_nat_pdb_dist12;

int* ibead_lj_non_nat;
int* jbead_lj_non_nat;
int* itype_lj_non_nat;
int* jtype_lj_non_nat;

// neighbor / cell list
int* ibead_neighbor_list_att;
int* jbead_neighbor_list_att;
int* itype_neighbor_list_att;
int* jtype_neighbor_list_att;

double* nl_lj_nat_pdb_dist;
double* nl_lj_nat_pdb_dist2;
double* nl_lj_nat_pdb_dist6;
double* nl_lj_nat_pdb_dist12;

int* ibead_neighbor_list_rep;
int* jbead_neighbor_list_rep;
int* itype_neighbor_list_rep;
int* jtype_neighbor_list_rep;

// pair list
int* ibead_pair_list_att;
int* jbead_pair_list_att;
int* itype_pair_list_att;
int* jtype_pair_list_att;

double* pl_lj_nat_pdb_dist;
double* pl_lj_nat_pdb_dist2;
double* pl_lj_nat_pdb_dist6;
double* pl_lj_nat_pdb_dist12;

int* ibead_pair_list_rep;
int* jbead_pair_list_rep;
int* itype_pair_list_rep;
int* jtype_pair_list_rep;

int lj_rna_rna_allocated = 0;

// coordinates and associated params

int nbead;
coord* pos;
coord* unc_pos; // uncorrected positions
coord* vel;
coord* force;
int pos_allocated = 0;
int vel_allocated = 0;
int force_allocated = 0;
int unc_pos_allocated = 0;

// native info

int* rna_base; // array which indicates whether or not a bead is a base
int rna_base_allocated;
int* rna_phosphate;
int rna_phosphate_allocated;

// miscellaneous run paramaters;

Ran_Gen generator; // random number generator
int run;
int restart = 0; // default is to start a new simulation
int rgen_restart = 0; // default don't restart random number generator
int sim_type = 1; // integration scheme; default is underdamped
double T; // temperature
int neighborlist = 0; // neighbor list cutoff method?
int celllist = 0; // cell list cutoff method?
double boxl; // Length of an edge of the simulation box
double ncell;
double lcell;
double zeta; // friction coefficient
double nstep; // number of steps to take
double istep_restart = 0.0;
int nup;
int inlup;
int nnlup;
double h; // time step
double halfh;
double a1; // a1,a2,a3,a4,a5 are used for integration
double a2;
double a3;
double a4;
double a5;
char ufname[mwdsize+1];
char rcfname[mwdsize+1];
char cfname[mwdsize+1];
char unccfname[mwdsize+1];
char vfname[mwdsize+1];
char binfname[mwdsize+1];
char uncbinfname[mwdsize+1];
char iccnfigfname[mwdsize+1];
int binsave = 1; // default will save trajectory
const double pi = acos(-1);

// force and pot stuff

int nforce_term = 4; // ran,bnds,angs,vdw -- default is that tension is off
int force_term_on[mforce_term+1] = { 0, 1, 1, 1, 0,
				     0, 1, 0, 0, 0, 0 };
force_term_Ptr force_term[mforce_term+1];

int npot_term = 3; // bnds,angs,vdw
int pot_term_on[mpot_term+1] = { 0, 1, 1, 0, 0,
				 1, 0, 0, 0, 0, 0 };
pot_term_Ptr pot_term[mpot_term+1];

//observables

double e_bnd,e_ang,e_tor,e_stack,e_elec,e_ang_ss;
double e_vdw_rr,e_vdw_rr_att,e_vdw_rr_rep;
double e_vdw_cc,e_vdw_rc,e_vdw_rc_rep,e_vdw_rc_att;
double rna_etot,system_etot;
double chi;
double Q;
int contct_nat;
int contct_tot;
double end2endsq;
double rgsq;
double kinT;

int main(int argc,char* argv[])
{

  using namespace std;

  if( argc<2 ) {
    cerr << "Usage: " << argv[0] <<  " < input_file >" << endl;
    exit(-1);    
  }
  time_t tm0 = time(0); // wall time at this point
  cout << "CURRENT TIME IS: " << ctime(&tm0);
  if( getcwd(pathname,MAXPATHLEN)==NULL ) {
    cerr << "PROBLEM GETTING PATH" << endl;
  } else {
    cout << "CURRENT WORKING DIRECTORY: " << pathname << endl;
  }

  alloc_arrays(); // allocates certain arrays and initializes some variables
  read_input(argv[1]); // read input file

  clock_t ck0 = clock(); // clock ticks at this point
  ex_cmds(); // perform commands (simulation)

  // time stats
  time_t tm1 = time(0);
  clock_t ck1 = clock();
  cout << "+-------------------+" << endl;
  cout << "| Simulation Stats: |" << endl;
  cout << "+-------------------+" << endl;
  cout << "Wall Time              : " << difftime(tm1,tm0) << " sec" << endl;
  cout << "Total Computation Time : " << float(ck1-ck0)/CLOCKS_PER_SEC << " sec" << endl;
  cout << "Computation Rate       : " << float(ck1-ck0)/CLOCKS_PER_SEC/nstep << " sec / timestep" << endl;
  cout << "CURRENT TIME IS        : " << ctime(&tm1);

  return 0;

}

void alloc_arrays()
{

  using namespace std;

  // bonds

  k_bnd = 20.0;
  R0 = 2.0; // = 0.4*a
  R0sq = R0*R0;
  e_bnd_coeff = k_bnd*R0sq/2.0; // SOP model
  nbnd = 1529;
  ibead_bnd = new int[nbnd+1];
  jbead_bnd = new int[nbnd+1];
  pdb_dist = new double[nbnd+1];
  bnds_allocated = 1;

  // angles

  k_ang = 20.0;
  e_ang_coeff = k_ang/2.0;
  nang = 1528;
  ibead_ang = new int[nang+1];
  jbead_ang = new int[nang+1];
  kbead_ang = new int[nang+1];
  pdb_ang = new double[nang+1];
  angs_allocated = 1;

  sigma_ss = 3.5; // = 0.76*a
  sigma_ss6 = pow(sigma_ss,6.0);
  epsilon_ss = 1.0;
  e_ang_ss_coeff = epsilon_ss*sigma_ss6;
  f_ang_ss_coeff = 6.0*e_ang_ss_coeff;

  // rna-rna vdw

  ncon_att = 8996;
  ncon_rep = 1157632;
  // neighbor list
  nnl_att = 0;
  nnl_rep = 0;
  // pair list
  nil_att = 0;
  nil_rep = 0;

  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      sigma_rep2[i][j] = sigma_rep[i][j] * sigma_rep[i][j];
      sigma_rep6[i][j] = sigma_rep2[i][j] * sigma_rep2[i][j] * sigma_rep2[i][j];
      sigma_rep12[i][j] = sigma_rep6[i][j] * sigma_rep6[i][j];
    }
  }

  ibead_lj_nat = new int[ncon_att+1];
  jbead_lj_nat = new int[ncon_att+1];
  itype_lj_nat = new int[ncon_att+1];
  jtype_lj_nat = new int[ncon_att+1];
  lj_nat_pdb_dist = new double[ncon_att+1];
  lj_nat_pdb_dist2 = new double[ncon_att+1];
  lj_nat_pdb_dist6 = new double[ncon_att+1];
  lj_nat_pdb_dist12 = new double[ncon_att+1];
  ibead_lj_non_nat = new int[ncon_rep+1];
  jbead_lj_non_nat = new int[ncon_rep+1];
  itype_lj_non_nat = new int[ncon_rep+1];
  jtype_lj_non_nat = new int[ncon_rep+1];

  ibead_neighbor_list_att = new int[ncon_att+1];
  jbead_neighbor_list_att = new int[ncon_att+1];
  itype_neighbor_list_att = new int[ncon_att+1];
  jtype_neighbor_list_att = new int[ncon_att+1];
  nl_lj_nat_pdb_dist = new double[ncon_att+1];
  nl_lj_nat_pdb_dist2 = new double[ncon_att+1];
  nl_lj_nat_pdb_dist6 = new double[ncon_att+1];
  nl_lj_nat_pdb_dist12 = new double[ncon_att+1];
  ibead_neighbor_list_rep = new int[ncon_rep+1];
  jbead_neighbor_list_rep = new int[ncon_rep+1];
  itype_neighbor_list_rep = new int[ncon_rep+1];
  jtype_neighbor_list_rep = new int[ncon_rep+1];

  ibead_pair_list_att = new int[ncon_att+1];
  jbead_pair_list_att = new int[ncon_att+1];
  itype_pair_list_att = new int[ncon_att+1];
  jtype_pair_list_att = new int[ncon_att+1];
  pl_lj_nat_pdb_dist = new double[ncon_att+1];
  pl_lj_nat_pdb_dist2 = new double[ncon_att+1];
  pl_lj_nat_pdb_dist6 = new double[ncon_att+1];
  pl_lj_nat_pdb_dist12 = new double[ncon_att+1];
  ibead_pair_list_rep = new int[ncon_rep+1];
  jbead_pair_list_rep = new int[ncon_rep+1];
  itype_pair_list_rep = new int[ncon_rep+1];
  jtype_pair_list_rep = new int[ncon_rep+1];

  lj_rna_rna_allocated = 1;

  // coordinates
  
  nbead = 1530;
  pos = new coord[nbead+1];
  unc_pos = new coord[nbead+1];
  vel = new coord[nbead+1];
  force = new coord[nbead+1];
  rna_base = new int [nbead+1];
  rna_phosphate = new int [nbead+1];
  pos_allocated = 1;
  unc_pos_allocated = 1;
  vel_allocated = 1;
  force_allocated = 1;
  rna_base_allocated = 1;
  rna_phosphate_allocated = 1;

  // miscellaneous run parameters

  run = 1;
  generator.set_seed(-100-run);
  T = 0.6; // kcal/mol

  neighborlist = 0; // neighbor list cutoff method?
  celllist = 0; // cell list cutoff method?
  boxl = 500.0;
  ncell = 55.0;
  lcell = boxl / ncell;
  zeta = 5.0e-2; // 0.05*tau^{-1} = friction coeff
  nstep = 5e7;
  nup = 1000;
  nnlup = 50; // neighbor list update frequency
  h = 2.5e-3;
  halfh = h/2.0;
  a1 = h*(1.0-zeta*halfh);
  a2 = h*halfh;
  a3 = (1.0-h*zeta/2.0+(h*zeta)*(h*zeta)/4.0)/h;
  a4 = halfh*(1.0-h*zeta/2.0);
  a5 = h/zeta;
  strcpy(ufname,"update.out");
  strcpy(rcfname,"restart_c.dat");
  strcpy(cfname,"coord.out");
  strcpy(unccfname,"unccoord.out");
  strcpy(vfname,"veloc.out");
  strcpy(binfname,"traj.bin");
  strcpy(uncbinfname,"traj_uncorrected.bin");

}

void read_input(const char* const ifile)
{

  using namespace std;

  ifstream in;
  char line[1024];
  char* tokPtr;
  char term = ';'; // terminates a command
  int newcmd = 1;
  int icmd;
  int nopt_tot = 0;
  int iopt;

  ncmd = 0;
  in.clear();
  in.open(ifile,ios::in);
  while(1) {
    in.getline(line,1024);
    if( in.eof() ) break;
    tokPtr = strtok(line," ");
    if( strchr(tokPtr,term)!=NULL ) { ncmd++; }
    while( tokPtr = strtok(NULL," ") ) {
      if( strchr(tokPtr,term)!=NULL ) { ncmd++; }
    }      
  }
  in.close();

  //  cout << "NUMBER OF COMMANDS: " << ncmd << endl;

  in.clear();
  in.open(ifile,ios::in);
  icmd = 0;
  while(1) {
    in.getline(line,1024);
    if( in.eof() ) break;
    tokPtr = strtok(line," ");;
    if( newcmd ) { 
      icmd++; 
      strcpy( cmd[icmd],tokPtr ); 
      opt_ptr[icmd] = nopt_tot+1; 
      newcmd = 0;
    } else { 
      nopt_tot++; 
      strcpy( opt[nopt_tot],tokPtr ); 
    }
    if( strchr(tokPtr,term)!=NULL ) {
      newcmd = 1;
    }
    while( tokPtr = strtok(NULL," ") ) {
      if( newcmd ) { 
        icmd++; 
        strcpy( cmd[icmd],tokPtr ); 
        opt_ptr[icmd] = nopt_tot+1; 
        newcmd = 0;
      } else { 
        nopt_tot++; 
        strcpy( opt[nopt_tot],tokPtr ); 
      }
      if( strchr(tokPtr,term)!=NULL ) {
        newcmd = 1;
      }
    }
  }
  opt_ptr[ncmd+1] = nopt_tot + 1;
  in.close();

  for( int icmd=1; icmd<=ncmd; icmd++ ) {
    for( int i=0; i<strlen(cmd[icmd]); i++ ) {
      if( cmd[icmd][i] == term ) cmd[icmd][i] = '\0';
    }
    //    cout << "COMMAND[" << icmd << "]: " << cmd[icmd] << endl;
    for( int iopt = opt_ptr[icmd]; iopt < opt_ptr[icmd+1]; iopt++ ) {
      for( int i=0; i<strlen(opt[iopt]); i++ ) {
        if( opt[iopt][i] == ';' ) opt[iopt][i] = '\0';
      }
      //      cout << opt[iopt] << endl;
    }
  }
  
}

void ex_cmds()
{

  using namespace std;

  char oline[1024];
  int iopt;

  for( int i=1; i<=ncmd; i++ ) {
     // read data
     if( !strcmp(cmd[i],"load") ) { load(i); }
     // set parameters
     else if( !strcmp(cmd[i],"set") ) { set_params(i); }
     // run simulation
     else if( !strcmp(cmd[i],"run") ) { simulation_ctrl(); }
     // ???
     else {};
  }

}

void set_params(int icmd)
{

  using namespace std;
  char oline[1024];
  int iopt;

  if( !strcmp(opt[opt_ptr[icmd]],"dynamics") ) { // set the type of simulation
    if( !strcmp(opt[opt_ptr[icmd]+1],"underdamped") ) {
      sim_type = 1; // low-friction limit for sampling
      h = 2.5e-3;
      halfh = h/2.0;
      a1 = h*(1.0-zeta*halfh);
      a2 = h*halfh;
      a3 = (1.0-h*zeta/2.0+(h*zeta)*(h*zeta)/4.0)/h;
      a4 = halfh*(1.0-h*zeta/2.0);
    } else if( !strcmp(opt[opt_ptr[icmd]+1],"overdamped") ) {
      sim_type = 2; // hi-friction limit for kinetics
      h = 0.02;
      a5 = h/zeta;
    }
  } else if( !strcmp(opt[opt_ptr[icmd]],"temp") ) { // set the temperature
    set_temp(atof(opt[opt_ptr[icmd]+1]));

  } else if( !strcmp(opt[opt_ptr[icmd]],"nstep") ) { // # of steps
    nstep = atof(opt[opt_ptr[icmd]+1]);
  } else if( !strcmp(opt[opt_ptr[icmd]],"istep_restart") ) { // where to restart from
    istep_restart = atof(opt[opt_ptr[icmd]+1]);
  } else if( !strcmp(opt[opt_ptr[icmd]],"nup") ) { // # of steps before an update
    nup = atoi(opt[opt_ptr[icmd]+1]);

  } else if( !strcmp(opt[opt_ptr[icmd]],"run") ) { // set current run
    run = atoi((opt[opt_ptr[icmd]+1]));
    generator.set_seed(-100-run);

  } else if( !strcmp(opt[opt_ptr[icmd]],"ufname") ) { // set update file name
    strcpy(ufname,opt[opt_ptr[icmd]+1]);
  } else if( !strcmp(opt[opt_ptr[icmd]],"rcfname") ) { // set restart coordinate file name
    strcpy(rcfname,opt[opt_ptr[icmd]+1]);
  } else if( !strcmp(opt[opt_ptr[icmd]],"cfname") ) { // set save coordinate file name
    strcpy(cfname,opt[opt_ptr[icmd]+1]);
  } else if( !strcmp(opt[opt_ptr[icmd]],"rgenfname") ) { // set random generator file name
    generator.set_fname(opt[opt_ptr[icmd]+1]);
  } else if( !strcmp(opt[opt_ptr[icmd]],"unccfname") ) { // set save coordinate file name
    strcpy(unccfname,opt[opt_ptr[icmd]+1]);
  } else if( !strcmp(opt[opt_ptr[icmd]],"vfname") ) { // set save velocity file name
    strcpy(vfname,opt[opt_ptr[icmd]+1]);
  } else if( !strcmp(opt[opt_ptr[icmd]],"binfname") ) { // set save trajectory file name
    strcpy(binfname,opt[opt_ptr[icmd]+1]);
  } else if( !strcmp(opt[opt_ptr[icmd]],"uncbinfname") ) { // set save trajectory file name
    strcpy(uncbinfname,opt[opt_ptr[icmd]+1]);

  } else if( !strcmp(opt[opt_ptr[icmd]],"cutofftype") ) { // neighbor list on or off?
    if( !strcmp(opt[opt_ptr[icmd]+1],"neighborlist" ) ) { neighborlist = 1; }
    else if( !strcmp(opt[opt_ptr[icmd]+1],"celllist" ) ) { celllist = 1; }
    else { }

  } else if( !strcmp(opt[opt_ptr[icmd]],"nnlup") ) { // neighbor / cell list update frequency
    nnlup = atoi(opt[opt_ptr[icmd]+1]);

  } else if( !strcmp(opt[opt_ptr[icmd]],"boxl") ) { // box length for pbc
    boxl = atof(opt[opt_ptr[icmd]+1]);
  } else if( !strcmp(opt[opt_ptr[icmd]],"ncell") ) { // number of cells along box length
    ncell = atof(opt[opt_ptr[icmd]+1]);
    lcell = boxl / ncell;

  } else if( !strcmp(opt[opt_ptr[icmd]],"restart") ) { // restart on or off?
    if( !strcmp(opt[opt_ptr[icmd]+1],"on" ) ) { restart = 1; }
    else { restart = 0; }
  } else if( !strcmp(opt[opt_ptr[icmd]],"rgen_restart") ) { // restart the generator?
    if( !strcmp(opt[opt_ptr[icmd]+1],"on" ) ) { rgen_restart = 1; }
    else { rgen_restart = 0; }
  } else if( !strcmp(opt[opt_ptr[icmd]],"t_step") ) {
    h = atof((opt[opt_ptr[icmd]+1]));
    halfh = h/2.0;
    a1 = h*(1.0-zeta*halfh);
    a2 = h*halfh;
    a3 = (1.0-h*zeta/2.0+(h*zeta)*(h*zeta)/4.0)/h;
    a4 = halfh*(1.0-h*zeta/2.0);
    a5 = h/zeta;
  } else if( !strcmp(opt[opt_ptr[icmd]],"zeta") ) { // friction coefficient
    if (sim_type == 1) h = 2.5e-3; else if (sim_type == 2) h = 0.02;
    zeta = atof((opt[opt_ptr[icmd]+1]));
    a1 = h*(1.0-zeta*halfh);
    a3 = (1.0-h*zeta/2.0+(h*zeta)*(h*zeta)/4.0)/h;
    a4 = halfh*(1.0-h*zeta/2.0);
    a5 = h/zeta;
  } else {};
  
}

void set_temp(double temp)
{
  using namespace std;

  T = temp;
}

void simulation_ctrl()
{

  using namespace std;

  switch( sim_type ) {
  case 1:
    underdamped_ctrl();
    break;
  case 2:
    overdamped_ctrl();
    break;
  default:
    cerr << "UNRECOGNIZED SIM_TYPE!" << endl;
    exit(-1);
  }

}

void underdamped_ctrl()
{

  using namespace std;

  char oline[2048];
  double istep = 1.0;
  int iup = 1;
  int inlup = 1;
  ofstream out(ufname,ios::out|ios::app);
  static int first_time = 1;

  coord* incr = new coord[nbead+1];

  if( (!restart)&&first_time ) { // zero out the velocities and forces
    for( int i=1; i<=nbead; i++ ) {
      vel[i].x = 0.0;
      vel[i].y = 0.0;
      vel[i].z = 0.0;
      force[i].x = 0.0;
      force[i].y = 0.0;
      force[i].z = 0.0;
    }
  }

  print_sim_params();

  if (neighborlist == 1) {
    update_neighbor_list();
    update_pair_list();
  } else if (celllist == 1) {
    update_cell_list();
    update_pair_list();
  }

  set_potential();
  set_forces();
  
  char line[2048];

  if( restart ) {
    load_coords(cfname,unccfname);
    load_vels(vfname);
    istep = istep_restart + 1.0;
  }
  
  if( rgen_restart ) {
    generator.restart();
  }

  if( first_time ) {

    energy_eval();
    force_eval();

  }

  if( binsave ) {
    if( (first_time)&&(!rgen_restart) ) {
      record_traj(binfname,uncbinfname);     
    }
    while( istep <= nstep ) {
      
      // compute pair separation list
      if ((inlup % nnlup) == 0) {
        if (neighborlist == 1) {
          update_neighbor_list();
        } else if (celllist == 1) {
          update_cell_list();
        }
	//	fprintf(stderr, "(%.0lf) neighbor list: (%d/%d)\n", istep, nnl_att, nnl_rep);
        inlup = 0;
      }
      inlup++;

      if (neighborlist == 1 || celllist == 1) {
        update_pair_list();
//	fprintf(stderr, "(%.0lf) pair list: (%d/%d)\n", istep, nil_att, nil_rep);
      }

      underdamped_iteration(incr);
      if( !(iup%nup) ) { // updates
	energy_eval();
	calculate_observables(incr);
        sprintf(oline,"%.0lf %f %f %f %f %f %f %f %d %f",
                istep,T,kinT,e_bnd,e_ang_ss,e_vdw_rr,rna_etot,
                Q,contct_nat,rgsq);
	out << oline << endl;
	iup = 0;
	record_traj(binfname,uncbinfname);     
	save_coords(cfname,unccfname);
	save_vels(vfname);
	generator.save_state();
      }
      istep += 1.0;
      iup++;
      
    }
    out.close();
  }
  
  if( first_time ) first_time = 0;
  
  delete [] incr;

  return;

}

void print_sim_params() {

  using namespace std;

  char oline[2048];

  cout << endl;
  sprintf(oline,"+------------------------+");
  cout << oline << endl;
  sprintf(oline,"| Simulation Parameters: |");
  cout << oline << endl;
  sprintf(oline,"+------------------------+");
  cout << oline << endl;

  if (sim_type == 1) {
    sprintf(oline,"Simulation Type                   : %s", "Underdamped");
    cout << oline << endl;
  } else if (sim_type == 2) {
    sprintf(oline,"Simulation Type                   : %s", "Overdamped");
    cout << oline << endl;
  } else {
    cerr << "UNRECOGNIZED SIMULATION TYPE!" << endl;
    exit(-1);
  }

  sprintf(oline,"Simulation Temperature            : %.3f",T);
  cout << oline << endl;

  sprintf(oline,"Start Time Step                   : %.0lf", istep_restart);
  cout << oline << endl;

  sprintf(oline,"Final Time Step                   : %.0lf", nstep);
  cout << oline << endl;

  sprintf(oline,"Output Frequency                  : %d", nup);
  cout << oline << endl;

  sprintf(oline,"Friction Coefficient              : %.0e", zeta);
  cout << oline << endl;

  sprintf(oline,"PBC Box Length                    : %.1f", boxl);
  cout << oline << endl;

  if (neighborlist == 1) {
    sprintf(oline,"Long-range Cutoff Type            : %s", "Neighbor List");
    cout << oline << endl;
    sprintf(oline,"Neighbor List Update Frequency    : %d", nnlup);
    cout << oline << endl;
  } else if (celllist == 1) {
    sprintf(oline,"Long-range Cutoff Type            : %s", "Cell List");
    cout << oline << endl;
    sprintf(oline,"Cell List Update Frequency        : %d", nnlup);
    cout << oline << endl;

    sprintf(oline,"Number of Cells Each Dimension    : %.0lf", ncell);
    cout << oline << endl;
  } else {
    sprintf(oline,"Long-range Cutoff Type            : %s", "None");
    cout << oline << endl;
  }

  cout << endl;

}

void save_coords(char* fname,char* fname2)
{

  using namespace std;

  char oline[1024];
  ofstream ofile;
  ofstream ofile2;

  ofile.open(fname,ios::out);
  ofile2.open(fname2,ios::out);
  for( int i=1; i<=nbead; i++ ) {
    sprintf(oline,"%d %f %f %f",i,pos[i].x,
	    pos[i].y,pos[i].z);
    ofile << oline << endl;
    sprintf(oline,"%d %f %f %f",i,unc_pos[i].x,
	    unc_pos[i].y,unc_pos[i].z);
    ofile2 << oline << endl;
  }
  ofile.close();
  ofile2.close();
  
}

void load_coords(char* fname,char* fname2)
{

  using namespace std;

  char iline[1024];
  ifstream ifile;
  ifstream ifile2;
  char* tokPtr;

  ifile.clear();
  ifile2.clear();
  ifile.open(fname,ios::in);
  ifile2.open(fname2,ios::in);
  for( int i=1; i<=nbead; i++ ) {
    ifile.getline(iline,1024);
    tokPtr = strtok(iline," ");
    tokPtr = strtok(NULL," ");
    pos[i].x = atof(tokPtr);
    tokPtr = strtok(NULL," ");
    pos[i].y = atof(tokPtr);
    tokPtr = strtok(NULL," ");
    pos[i].z = atof(tokPtr);

    ifile2.getline(iline,1024);
    tokPtr = strtok(iline," ");
    tokPtr = strtok(NULL," ");
    unc_pos[i].x = atof(tokPtr);
    tokPtr = strtok(NULL," ");
    unc_pos[i].y = atof(tokPtr);
    tokPtr = strtok(NULL," ");
    unc_pos[i].z = atof(tokPtr);
  }
  ifile.close();
  ifile2.close();

}

void load_vels(char* fname)
{

  using namespace std;

  char iline[1024];
  ifstream ifile;
  char* tokPtr;

  ifile.clear();
  ifile.open(fname,ios::in);
  for( int i=1; i<=nbead; i++ ) {
    ifile.getline(iline,1024);
    tokPtr = strtok(iline," ");
    tokPtr = strtok(NULL," ");
    vel[i].x = atof(tokPtr);
    tokPtr = strtok(NULL," ");
    vel[i].y = atof(tokPtr);
    tokPtr = strtok(NULL," ");
    vel[i].z = atof(tokPtr);

  }
  ifile.close();

}

void save_unccoords(char* fname)
{

  using namespace std;

  char oline[1024];
  ofstream ofile;

  ofile.open(fname,ios::out);
  for( int i=1; i<=nbead; i++ ) {
    sprintf(oline,"%d %f %f %f",i,unc_pos[i].x,
	    unc_pos[i].y,unc_pos[i].z);
    ofile << oline << endl;
  }
  ofile.close();

}

void save_vels(char* fname)
{

  using namespace std;

  char oline[1024];
  ofstream ofile;

  ofile.open(fname,ios::out);
  for( int i=1; i<=nbead; i++ ) {
    sprintf(oline,"%d %f %f %f",i,vel[i].x,
	    vel[i].y,vel[i].z);
    ofile << oline << endl;
  }
  ofile.close();

}

void record_traj(char* fname,char* fname2)
{

  using namespace std;

  char oline[1024];
  char oline2[1024];
  ofstream trajfile;
  ofstream trajfile2;

  trajfile.open(fname,ios::out | ios::app);
  trajfile2.open(fname2,ios::out | ios::app);

  for( int i=1; i<=nbead; i++ ) {
    sprintf(oline,"%f %f %f",pos[i].x,pos[i].y,pos[i].z);
    sprintf(oline2,"%f %f %f",unc_pos[i].x,unc_pos[i].y,unc_pos[i].z);
    trajfile << oline << endl;
    trajfile2 << oline2 << endl;

  }

  trajfile.close();
  trajfile2.close();

} 

void calculate_observables(coord* increment)
{

  using namespace std;

  char oline[1024];
  double dx,dy,dz,d;
  static const double tol = 1.0; // tolerance for chi distances
  static const double chinorm = (double(nbead*nbead)-5.0*double(nbead)+6.0)/2.0;
  double sumvsq;
  int nchi;
  int ibead, jbead;
  int itype, jtype;
  float r_ij;
  char line[2048];

  // chi, contct_nat, contct_tot, Q
  
  contct_nat = 0;
  for( int i=1; i<=ncon_att; i++ ) {

    ibead = ibead_lj_nat[i];
    jbead = jbead_lj_nat[i];
    r_ij = lj_nat_pdb_dist[i];
    itype = itype_lj_nat[i];
    jtype = jtype_lj_nat[i];

    dx = unc_pos[ibead].x-unc_pos[jbead].x;
    dy = unc_pos[ibead].y-unc_pos[jbead].y;
    dz = unc_pos[ibead].z-unc_pos[jbead].z;

    dx -= boxl*rnd(dx/boxl);
    dy -= boxl*rnd(dy/boxl);
    dz -= boxl*rnd(dz/boxl);

    d = sqrt( dx*dx+dy*dy+dz*dz );
    if(d/r_ij < 1.25) {
      contct_nat++;
    }
  }
  Q = double(contct_nat)/ncon_att;

  
  // rgsq
  
  rgsq = 0.0;
  for( int i=1; i<=nbead-1; i++ ) {
    for( int j=i+1; j<=nbead; j++ ) {
      dx = unc_pos[i].x-unc_pos[j].x;
      dy = unc_pos[i].y-unc_pos[j].y;
      dz = unc_pos[i].z-unc_pos[j].z;
      dx -= boxl*rnd(dx/boxl);
      dy -= boxl*rnd(dy/boxl);
      dz -= boxl*rnd(dz/boxl);

      rgsq += (dx*dx+dy*dy+dz*dz);
    }
  }
  rgsq /= double(nbead*nbead);
  
  // kinT

  if( sim_type == 1 ) {
    sumvsq = 0.0;
    for( int i=1; i<=nbead; i++ ) {
      sumvsq += vel[i].x*vel[i].x
	+ vel[i].y*vel[i].y
	+ vel[i].z*vel[i].z;
    }
    kinT = sumvsq/(3.0*double(nbead));
  } else if( sim_type == 2 ) {
    sumvsq = 0.0;
    for( int i=1; i<=nbead; i++ ) {
      sumvsq += increment[i].x*increment[i].x +
	increment[i].y*increment[i].y +
	increment[i].z*increment[i].z;
    }
    sumvsq *= zeta/(2.0*h);
    kinT = sumvsq/(3.0*double(nbead));
  } else {}
  
  
}

void underdamped_iteration(coord* incr)
{
  
  using namespace std;
  
  static const double eps = 1.0e-5;

  for( int i=1; i<=nbead; i++ ) {
    
    // compute position increments
    
    incr[i].x = a1*vel[i].x + a2*force[i].x;
    incr[i].y = a1*vel[i].y + a2*force[i].y;
    incr[i].z = a1*vel[i].z + a2*force[i].z;

    // update bead positions

    pos[i].x += incr[i].x;
    pos[i].y += incr[i].y;
    pos[i].z += incr[i].z;

    pos[i].x -= boxl*rnd(pos[i].x/boxl);
    pos[i].y -= boxl*rnd(pos[i].y/boxl);
    pos[i].z -= boxl*rnd(pos[i].z/boxl);

    unc_pos[i].x += incr[i].x;
    unc_pos[i].y += incr[i].y;
    unc_pos[i].z += incr[i].z;

  }

  // force_update

  force_eval();

  if( T < eps ) return; // don't update velocities for steepest descent

  // update_velocities

  for( int i=1; i<=nbead; i++ ) {

    // compute velocity increments
    
    vel[i].x = a3*incr[i].x + a4*force[i].x;
    vel[i].y = a3*incr[i].y + a4*force[i].y;
    vel[i].z = a3*incr[i].z + a4*force[i].z;
    
  }
  
}

void overdamped_iteration(coord* incr)
{
   using namespace std;

   for( int i=1; i<=nbead; i++ ) {

      // compute position increments

      incr[i].x = a5*force[i].x;
      incr[i].y = a5*force[i].y;
      incr[i].z = a5*force[i].z;

      // update bead positions

      unc_pos[i].x += incr[i].x;
      unc_pos[i].y += incr[i].y;
      unc_pos[i].z += incr[i].z;

      pos[i].x += incr[i].x;
      pos[i].y += incr[i].y;
      pos[i].z += incr[i].z;

      pos[i].x -= boxl*rnd(pos[i].x/boxl);
      pos[i].y -= boxl*rnd(pos[i].y/boxl);
      pos[i].z -= boxl*rnd(pos[i].z/boxl);

   }

   // force_update

   force_eval();

}

void energy_eval()
{
  
  using namespace std;
  char oline[1024];

  for( int i=1; i<=npot_term; i++ ) {
    pot_term[i]();
  }

  rna_etot = e_bnd + e_ang_ss + e_vdw_rr;
  system_etot = rna_etot + e_vdw_rc + e_vdw_cc;

}

void force_eval()
{

  using namespace std;
  char oline[1024];

  clear_forces();
  
  for( int i=1; i<=nforce_term; i++ ) {
    force_term[i]();
  }

}

void clear_forces() {
  
  using namespace std;

  for( int i=1; i<=nbead; i++ ) {
    force[i].x = 0.0;
    force[i].y = 0.0;
    force[i].z = 0.0;
  }
  
}

void set_potential() {

  using namespace std;
  
  int iterm;
  
  iterm = 0;
  for( int i=1; i<=mpot_term; i++ ) {
    switch(i) {
    case 1:
      if( pot_term_on[i] ) {
	pot_term[++iterm] = &fene_energy;
      }
      break;
    case 2:
      if( pot_term_on[i] ) {
	pot_term[++iterm] = &soft_sphere_angular_energy;
      }
      break;
    case 5:
      if( pot_term_on[i] ) {
	pot_term[++iterm] = &vdw_energy;
      }
      break;
    default:
      break;
    }
  }

}

void set_forces()
{
  
  using namespace std;

  int iterm;

  iterm = 0;
  for( int i=1; i<=mforce_term; i++ ) {
    switch(i) {
    case 1:
      if( force_term_on[i] ) {
	force_term[++iterm] = &random_force;
      }
      break;
    case 2:
      if( force_term_on[i] ) {
	force_term[++iterm] = &fene_forces;
      } 
      break;
    case 3:
      if( force_term_on[i] ) {
	force_term[++iterm] = &soft_sphere_angular_forces;
      }
      break;
    case 6:
      if( force_term_on[i] ) {
	force_term[++iterm] = &vdw_forces;
      }
      break;
    default:
      break;
    }
  }
  
}

void fene_energy()
{
  
  using namespace std;

  int ibead, jbead;
  double dx, dy, dz, d,dev;
  char line[2048];

  e_bnd = 0.0;
  for( int i=1; i<=nbnd; i++ ) {
    
    ibead = ibead_bnd[i];
    jbead = jbead_bnd[i];

    dx = unc_pos[jbead].x-unc_pos[ibead].x;
    dy = unc_pos[jbead].y-unc_pos[ibead].y;
    dz = unc_pos[jbead].z-unc_pos[ibead].z;


    // min images
    
    dx -= boxl*rnd(dx/boxl);
    dy -= boxl*rnd(dy/boxl);
    dz -= boxl*rnd(dz/boxl);

    d = sqrt(dx*dx+dy*dy+dz*dz);
    dev = d-pdb_dist[i];

    e_bnd += log1p(-dev*dev/R0sq); // log1p(x) = log(1-x)

  }
  
  e_bnd *= -e_bnd_coeff;
  
  return;
  
}

void soft_sphere_angular_energy()
{

  using namespace std;

  e_ang_ss = 0.0;
  int ibead, kbead;
  coord r_ik;
  double d,d6;

    for( int i=1; i<=nang; i++ ) {
      
      ibead = ibead_ang[i];
      kbead = kbead_ang[i];
    
      r_ik.x = unc_pos[kbead].x - unc_pos[ibead].x;
      r_ik.y = unc_pos[kbead].y - unc_pos[ibead].y;
      r_ik.z = unc_pos[kbead].z - unc_pos[ibead].z;

      // min images
      
      r_ik.x -= boxl*rnd(r_ik.x/boxl);
      r_ik.y -= boxl*rnd(r_ik.y/boxl);
      r_ik.z -= boxl*rnd(r_ik.z/boxl);

      d = sqrt(r_ik.x*r_ik.x + r_ik.y*r_ik.y + r_ik.z*r_ik.z);
      d6 = pow(d,6.0);

      e_ang_ss += e_ang_ss_coeff/d6;
  }
  
  return;
  
}

void vdw_energy()
{

  using namespace std;

  int ibead,jbead;
  int itype,jtype;
  double dx,dy,dz,d,d2,d6,d12;
  char line[2048];

  e_vdw_rr = 0.0;
  e_vdw_rr_att = 0.0;
  e_vdw_rr_rep = 0.0;
  e_vdw_cc = 0.0;
  e_vdw_rc = 0.0;

  for( int i=1; i<=nil_att; i++ ) {

    ibead = ibead_pair_list_att[i];
    jbead = jbead_pair_list_att[i];
    itype = itype_pair_list_att[i];
    jtype = jtype_pair_list_att[i];

    dx = unc_pos[jbead].x - unc_pos[ibead].x;
    dy = unc_pos[jbead].y - unc_pos[ibead].y;
    dz = unc_pos[jbead].z - unc_pos[ibead].z;

    // min images

    dx -= boxl*rnd(dx/boxl);
    dy -= boxl*rnd(dy/boxl);
    dz -= boxl*rnd(dz/boxl);
    
    d2 = dx*dx+dy*dy+dz*dz;
    d6 = d2*d2*d2;
    d12 = d6*d6;

    e_vdw_rr_att += coeff_att[itype][jtype] * (pl_lj_nat_pdb_dist12[i]/d12)-2.0*(pl_lj_nat_pdb_dist6[i]/d6);

  }

  for( int i=1; i<=nil_rep; i++ ) {
    
    ibead = ibead_pair_list_rep[i];
    jbead = jbead_pair_list_rep[i];
    itype = itype_pair_list_rep[i];
    jtype = jtype_pair_list_rep[i];

    dx = unc_pos[jbead].x - unc_pos[ibead].x;
    dy = unc_pos[jbead].y - unc_pos[ibead].y;
    dz = unc_pos[jbead].z - unc_pos[ibead].z;

    // min images

    dx -= boxl*rnd(dx/boxl);
    dy -= boxl*rnd(dy/boxl);
    dz -= boxl*rnd(dz/boxl);

    d2 = dx*dx+dy*dy+dz*dz;
    d6 = d2*d2*d2;
    d12 = d6*d6;    

    e_vdw_rr_rep += coeff_rep[itype][jtype] * (sigma_rep12[itype][jtype]/d12+sigma_rep6[itype][jtype]/d6);

  }

  e_vdw_rr = e_vdw_rr_att + e_vdw_rr_rep;

  return;

}

void overdamped_ctrl()
{

  using namespace std;

  char oline[2048];
  double istep = 1.0;
  int iup = 1;
  ofstream out(ufname,ios::out|ios::app);
  static int first_time = 1;

  coord* incr = new coord[nbead+1];

  if( (!restart)&&first_time ) { // zero out the velocities and forces
    for( int i=1; i<=nbead; i++ ) {
      vel[i].x = 0.0;
      vel[i].y = 0.0;
      vel[i].z = 0.0;
      force[i].x = 0.0;
      force[i].y = 0.0;
      force[i].z = 0.0;
    }
  }

  print_sim_params();

  if (neighborlist == 1) {
    update_neighbor_list();
    update_pair_list();
  } else if (celllist == 1) {
    update_cell_list();
    update_pair_list();
  }

  set_potential();
  set_forces();
  
  char line[2048];

  if( restart ) {
    load_coords(cfname,unccfname);
    //    load_vels(vfname);
    istep = istep_restart + 1.0;
  }
  
  if( rgen_restart ) {
    generator.restart();
  }

  if( first_time ) {

    energy_eval();
    force_eval();

  }

  if( binsave ) {
    if( (first_time)&&(!rgen_restart) ) {
      record_traj(binfname,uncbinfname);     
    }
    while( istep <= nstep ) {
      
      // compute pair separation list
      if ((inlup % nnlup) == 0) {
        if (neighborlist == 1) {
          update_neighbor_list();
        } else if (celllist == 1) {
          update_cell_list();
        }
	//	fprintf(stderr, "(%.0lf) neighbor list: (%d/%d)\n", istep, nnl_att, nnl_rep);
        inlup = 0;
      }
      inlup++;

      if (neighborlist == 1 || celllist == 1) {
        update_pair_list();
//	fprintf(stderr, "(%.0lf) pair list: (%d/%d)\n", istep, nil_att, nil_rep);
      }

      overdamped_iteration(incr);
      if( !(iup%nup) ) { // updates
	energy_eval();
	calculate_observables(incr);
        sprintf(oline,"%.0lf %f %f %f %f %f %f %f %d %f",
                istep,T,kinT,e_bnd,e_ang_ss,e_vdw_rr,rna_etot,
                Q,contct_nat,rgsq);
	out << oline << endl;
	iup = 0;
	record_traj(binfname,uncbinfname);     
	save_coords(cfname,unccfname);
	save_vels(vfname);
	generator.save_state();
      }
      istep += 1.0;
      iup++;
      
    }
    out.close();
  }
  
  if( first_time ) first_time = 0;
  
  delete [] incr;

  return;

}

coord::coord() {
  
  x = 0.0;
  y = 0.0;
  z = 0.0;

}

coord::~coord() {

  // nothing to do

}

void load(int icmd)
{
  
  using namespace std;

  ifstream in;
  char line[2048];
  char* tokPtr;
  int test;
  int test1,test2;
  int ncon_tot;
  int icon_att, icon_rep;
  int i,j,k,l;
  int ibead,jbead;
  int itype,jtype;
  double real_phi, ideal_phi;
  double r_ij;
  int istate;

  if( !strcmp(opt[opt_ptr[icmd]],"bonds") ) { // load bonds
    cout << "[Reading in bonds...]" << endl;
    in.clear();
    in.open(opt[opt_ptr[icmd]+1],ios::in); // open file
    in.getline(line,2048);
    tokPtr = strtok(line," ");
    tokPtr = strtok(NULL," ");
    nbnd = atoi(tokPtr); // read in number of bonds
    init_bonds(nbnd);
    for( int i=1; i<=nbnd; i++ ) {
      in.getline(line,2048);
      tokPtr = strtok(line," ");
      ibead_bnd[i] = atoi(tokPtr); // first bead index
      tokPtr = strtok(NULL," ");
      jbead_bnd[i] = atoi(tokPtr); // second bead index
      tokPtr = strtok(NULL," ");
      pdb_dist[i] = atof(tokPtr); // equilibrium distance (angstrom)
    }
    in.close(); // close file
    cout << "[Finished reading bonds (" << nbnd <<")]" << endl;
  } else if(!strcmp(opt[opt_ptr[icmd]],"angles")) { // load angles
    cout << "[Reading in angles...]" << endl;
    in.clear();
    in.open(opt[opt_ptr[icmd]+1],ios::in);
    in.getline(line,2048);
    tokPtr = strtok(line," ");
    tokPtr = strtok(NULL," ");
    nang = atoi(tokPtr); // read in number of angles
    init_angles(nang);
    for( int i=1; i<=nang; i++ ) {
      in.getline(line,2048);
      tokPtr = strtok(line," ");
      ibead_ang[i] = atoi(tokPtr); // first bead index
      tokPtr = strtok(NULL," ");
      jbead_ang[i] = atoi(tokPtr); // second bead index
      tokPtr = strtok(NULL," ");
      kbead_ang[i] = atoi(tokPtr); // third bead index
      tokPtr = strtok(NULL," ");
      pdb_ang[i] = atof(tokPtr); // equilibrium angle (radians) ; SOP -> dist between bead i,i+2
    }
    in.close();    
    cout << "[Finished reading angles (" << nang <<")]" << endl;
  } else if(!strcmp(opt[opt_ptr[icmd]],"vdw")) { // load rna-rna vdw
    cout << "[Reading in VDW interactions...]" << endl;
    in.clear();
    in.open(opt[opt_ptr[icmd]+1],ios::in);
    in.getline(line,2048);
    tokPtr = strtok(line," ");
    tokPtr = strtok(NULL," ");
    ncon_att = atoi(tokPtr);
    tokPtr = strtok(NULL," ");
    tokPtr = strtok(NULL," ");
    ncon_rep = atoi(tokPtr);
    init_lj(ncon_att,ncon_rep);
    ncon_tot = ncon_att + ncon_rep;
    icon_att = 0;
    icon_rep = 0;
    for( int i=1; i<=ncon_tot; i++ ) {
      in.getline(line,2048);
      tokPtr = strtok(line," ");
      ibead = atoi(tokPtr);
      tokPtr = strtok(NULL," ");
      jbead = atoi(tokPtr);
      tokPtr = strtok(NULL," ");
      r_ij = atof(tokPtr);
      tokPtr = strtok(NULL," ");
      itype = atoi(tokPtr);
      tokPtr = strtok(NULL," ");
      jtype = atoi(tokPtr);
      if (r_ij < rcut_nat[itype][jtype]) {
	icon_att++;
	ibead_lj_nat[icon_att] = ibead;
	jbead_lj_nat[icon_att] = jbead;
	itype_lj_nat[icon_att] = itype;
	jtype_lj_nat[icon_att] = jtype;
	lj_nat_pdb_dist[icon_att] = r_ij;
	lj_nat_pdb_dist2[icon_att] = r_ij*r_ij;
	lj_nat_pdb_dist6[icon_att] = lj_nat_pdb_dist2[icon_att]*
	  lj_nat_pdb_dist2[icon_att]*lj_nat_pdb_dist2[icon_att];
	lj_nat_pdb_dist12[icon_att] = lj_nat_pdb_dist6[icon_att]*
	  lj_nat_pdb_dist6[icon_att];

	nil_att++;
	ibead_pair_list_att[nil_att] = ibead;
	jbead_pair_list_att[nil_att] = jbead;
	itype_pair_list_att[nil_att] = itype;
	jtype_pair_list_att[nil_att] = jtype;
	pl_lj_nat_pdb_dist[nil_att] = r_ij;
	pl_lj_nat_pdb_dist2[nil_att] = lj_nat_pdb_dist2[icon_att];
	pl_lj_nat_pdb_dist6[nil_att] = lj_nat_pdb_dist6[icon_att];
	pl_lj_nat_pdb_dist12[nil_att] = lj_nat_pdb_dist12[icon_att];
      } else {
	icon_rep++;
	ibead_lj_non_nat[icon_rep] = ibead;
	jbead_lj_non_nat[icon_rep] = jbead;
	itype_lj_non_nat[icon_rep] = itype;
	jtype_lj_non_nat[icon_rep] = jtype;

	nil_rep++;
	ibead_pair_list_rep[nil_rep] = ibead;
	jbead_pair_list_rep[nil_rep] = jbead;
	itype_pair_list_rep[nil_rep] = itype;
	jtype_pair_list_rep[nil_rep] = jtype;
      }
    }
    in.close();
    cout << "[Finished reading VDW interactions (" << icon_att << "/" << icon_rep <<")]" << endl;
  } else if(!strcmp(opt[opt_ptr[icmd]],"init")) { // load init coordinates
    cout << "[Reading in initial coordinates...]" << endl;
    in.clear();
    in.open(opt[opt_ptr[icmd]+1],ios::in);
    in.getline(line,2048);
    tokPtr = strtok(line," ");
    tokPtr = strtok(NULL," ");
    nbead = atoi(tokPtr); // read in number of beads
    init_pos(nbead);
    for( int i=1; i<=nbead; i++ ) {
      in.getline(line,2048);
      tokPtr = strtok(line," ");
      tokPtr = strtok(NULL," ");
      pos[i].x = atof(tokPtr);
      unc_pos[i].x = pos[i].x;
      tokPtr = strtok(NULL," ");
      pos[i].y = atof(tokPtr);
      unc_pos[i].y = pos[i].y;
      tokPtr = strtok(NULL," ");
      pos[i].z = atof(tokPtr);
      unc_pos[i].z = pos[i].z;
    }
    in.close();
    cout << "[Finished reading initial coordinates (" << nbead << ")]" << endl;
  }
  
}

void release_bonds()
{

  using namespace std;

  delete [] ibead_bnd;
  delete [] jbead_bnd;
  delete [] pdb_dist;
  bnds_allocated = 0;

}

void init_bonds(int numbonds)
{
 
  using namespace std;

  nbnd = numbonds;
  ibead_bnd = new int[numbonds+1];
  jbead_bnd = new int[numbonds+1];
  pdb_dist = new double[numbonds+1];
  bnds_allocated = 1;

}

void release_angles()
{

  using namespace std;

  delete [] ibead_ang;
  delete [] jbead_ang;
  delete [] kbead_ang;
  delete [] pdb_ang;
  angs_allocated = 0;

}

void init_angles(int numangs)
{
 
  using namespace std;

  nang = numangs;
  ibead_ang = new int[numangs+1];
  jbead_ang = new int[numangs+1];
  kbead_ang = new int[numangs+1];
  pdb_ang = new double[numangs+1];
  angs_allocated = 1;

}

void release_lj()
{

  using namespace std;

  delete [] ibead_lj_nat;
  delete [] jbead_lj_nat;
  delete [] itype_lj_nat;
  delete [] jtype_lj_nat;
  delete [] lj_nat_pdb_dist;
  delete [] lj_nat_pdb_dist2;
  delete [] lj_nat_pdb_dist6;
  delete [] lj_nat_pdb_dist12;
  delete [] ibead_lj_non_nat;
  delete [] jbead_lj_non_nat;
  delete [] itype_lj_non_nat;
  delete [] jtype_lj_non_nat;

  delete [] ibead_neighbor_list_att;
  delete [] jbead_neighbor_list_att;
  delete [] itype_neighbor_list_att;
  delete [] jtype_neighbor_list_att;
  delete [] nl_lj_nat_pdb_dist;
  delete [] nl_lj_nat_pdb_dist2;
  delete [] nl_lj_nat_pdb_dist6;
  delete [] nl_lj_nat_pdb_dist12;
  delete [] ibead_neighbor_list_rep;
  delete [] jbead_neighbor_list_rep;
  delete [] itype_neighbor_list_rep;
  delete [] jtype_neighbor_list_rep;

  // pair list
  delete [] ibead_pair_list_att;
  delete [] jbead_pair_list_att;
  delete [] itype_pair_list_att;
  delete [] jtype_pair_list_att;
  delete [] pl_lj_nat_pdb_dist;
  delete [] pl_lj_nat_pdb_dist2;
  delete [] pl_lj_nat_pdb_dist6;
  delete [] pl_lj_nat_pdb_dist12;
  delete [] ibead_pair_list_rep;
  delete [] jbead_pair_list_rep;
  delete [] itype_pair_list_rep;
  delete [] jtype_pair_list_rep;

  lj_rna_rna_allocated = 0;

}

void init_lj(int numatt, int numrep ) 
{

  using namespace std;

  ncon_att = numatt;
  ncon_rep = numrep;
  ibead_lj_nat = new int[numatt+1];
  jbead_lj_nat = new int[numatt+1];
  itype_lj_nat = new int[numatt+1];
  jtype_lj_nat = new int[numatt+1];
  lj_nat_pdb_dist = new double[numatt+1];
  lj_nat_pdb_dist2 = new double[numatt+1];
  lj_nat_pdb_dist6 = new double[numatt+1];
  lj_nat_pdb_dist12 = new double[numatt+1];
  ibead_lj_non_nat = new int[numrep+1];
  jbead_lj_non_nat = new int[numrep+1];
  itype_lj_non_nat = new int[numrep+1];
  jtype_lj_non_nat = new int[numrep+1];

  ibead_neighbor_list_att = new int[numatt+1];
  jbead_neighbor_list_att = new int[numatt+1];
  itype_neighbor_list_att = new int[numatt+1];
  jtype_neighbor_list_att = new int[numatt+1];
  nl_lj_nat_pdb_dist = new double[numatt+1];
  nl_lj_nat_pdb_dist2 = new double[numatt+1];
  nl_lj_nat_pdb_dist6 = new double[numatt+1];
  nl_lj_nat_pdb_dist12 = new double[numatt+1];
  ibead_neighbor_list_rep = new int[numrep+1];
  jbead_neighbor_list_rep = new int[numrep+1];
  itype_neighbor_list_rep = new int[numrep+1];
  jtype_neighbor_list_rep = new int[numrep+1];

  ibead_pair_list_att = new int[numatt+1];
  jbead_pair_list_att = new int[numatt+1];
  itype_pair_list_att = new int[numatt+1];
  jtype_pair_list_att = new int[numatt+1];
  pl_lj_nat_pdb_dist = new double[numatt+1];
  pl_lj_nat_pdb_dist2 = new double[numatt+1];
  pl_lj_nat_pdb_dist6 = new double[numatt+1];
  pl_lj_nat_pdb_dist12 = new double[numatt+1];
  ibead_pair_list_rep = new int[numrep+1];
  jbead_pair_list_rep = new int[numrep+1];
  itype_pair_list_rep = new int[numrep+1];
  jtype_pair_list_rep = new int[numrep+1];

  lj_rna_rna_allocated = 1;

}

void init_pos(int nbead)
{

  using namespace std;

  unc_pos = new coord[nbead+1];
  pos = new coord[nbead+1];

  vel = new coord[nbead+1];
  force = new coord[nbead+1];

  pos_allocated = 1;
  unc_pos_allocated = 1;
  vel_allocated = 1;
  force_allocated = 1;
}

void release_pos()
{

  using namespace std;

  delete [] unc_pos;
  delete [] pos;

  delete [] vel;
  delete [] force;

  pos_allocated = 0;
  unc_pos_allocated = 0;
  vel_allocated = 0;
  force_allocated = 0;
}

void vdw_forces()
{

  using namespace std;

  char line[2048];

  int ibead,jbead;
  int itype,jtype;
  double dx,dy,dz,d,d2,d6,d12;
  double fx,fy,fz;
  double co1;
  const static double tol = 1.0e-7;
  double rep_tol;

  for( int i=1; i<=nil_att; i++ ) {

    ibead = ibead_pair_list_att[i];
    jbead = jbead_pair_list_att[i];
    itype = itype_pair_list_att[i];
    jtype = jtype_pair_list_att[i];

    dx = unc_pos[jbead].x - unc_pos[ibead].x;
    dy = unc_pos[jbead].y - unc_pos[ibead].y;
    dz = unc_pos[jbead].z - unc_pos[ibead].z;

    // min images

    dx -= boxl*rnd(dx/boxl);
    dy -= boxl*rnd(dy/boxl);
    dz -= boxl*rnd(dz/boxl);
    
    d2 = dx*dx+dy*dy+dz*dz;
    rep_tol = sigma_rep2[itype][jtype]*tol;
    if( d2 < tol*lj_nat_pdb_dist2[i] ) continue;
    d6 = d2*d2*d2;
    d12 = d6*d6;
    
    co1 = force_coeff_att[itype][jtype]/d2*((pl_lj_nat_pdb_dist12[i]/d12)-(pl_lj_nat_pdb_dist6[i]/d6));

    fx = co1*dx;
    fy = co1*dy;
    fz = co1*dz;

    force[ibead].x += fx;
    force[ibead].y += fy;
    force[ibead].z += fz;

    force[jbead].x -= fx;
    force[jbead].y -= fy;
    force[jbead].z -= fz;

  }

  for( int i=1; i<=nil_rep; i++ ) {

    ibead = ibead_pair_list_rep[i];
    jbead = jbead_pair_list_rep[i];
    itype = itype_pair_list_rep[i];
    jtype = jtype_pair_list_rep[i];

    dx = unc_pos[jbead].x - unc_pos[ibead].x;
    dy = unc_pos[jbead].y - unc_pos[ibead].y;
    dz = unc_pos[jbead].z - unc_pos[ibead].z;

    // min images

    dx -= boxl*rnd(dx/boxl);
    dy -= boxl*rnd(dy/boxl);
    dz -= boxl*rnd(dz/boxl);

    d2 = dx*dx+dy*dy+dz*dz;
    if( d2 <  rep_tol ) continue;    
    d6 = d2*d2*d2;
    d12 = d6*d6;    

    co1 = force_coeff_rep[itype][jtype]/d2*
      (2.0*sigma_rep12[itype][jtype]/d12+sigma_rep6[itype][jtype]/d6);
    
    fx = co1*dx;
    fy = co1*dy;
    fz = co1*dz;

    force[ibead].x += fx;
    force[ibead].y += fy;
    force[ibead].z += fz;

    force[jbead].x -= fx;
    force[jbead].y -= fy;
    force[jbead].z -= fz;

  }

}

void soft_sphere_angular_forces()
{

  using namespace std;

  char line[2048];

  int ibead,kbead;
  double dx,dy,dz,d,d8;
  double fx,fy,fz;
  double co1;

  for( int i=1; i<=nang; i++ ) {
    
      ibead = ibead_ang[i];
      kbead = kbead_ang[i];
    
      dx = unc_pos[kbead].x - unc_pos[ibead].x;
      dy = unc_pos[kbead].y - unc_pos[ibead].y;
      dz = unc_pos[kbead].z - unc_pos[ibead].z;
    
      // min images
    
      dx -= boxl*rnd(dx/boxl);
      dy -= boxl*rnd(dy/boxl);
      dz -= boxl*rnd(dz/boxl);
      
      d = sqrt(dx*dx+dy*dy+dz*dz);
      d8 = pow(d,8.0);
    
      co1 = f_ang_ss_coeff/d8;
      
      fx = co1*dx;
      fy = co1*dy;
      fz = co1*dz;
      
      force[ibead].x -= fx;
      force[ibead].y -= fy;
      force[ibead].z -= fz;
    
      force[kbead].x += fx;
      force[kbead].y += fy;
      force[kbead].z += fz;

  }
  
}

void fene_forces()
{

  using namespace std;


  int ibead, jbead;
  double dx, dy, dz, d, dev, dev2;
  double fx, fy, fz;
  double temp;

  char line[2048];

  for( int i=1; i<=nbnd; i++ ) {
    
    ibead = ibead_bnd[i];
    jbead = jbead_bnd[i];
    
    dx = unc_pos[jbead].x-unc_pos[ibead].x;
    dy = unc_pos[jbead].y-unc_pos[ibead].y;
    dz = unc_pos[jbead].z-unc_pos[ibead].z;

    dx -= boxl*rnd(dx/boxl);
    dy -= boxl*rnd(dy/boxl);
    dz -= boxl*rnd(dz/boxl);
    
    d = sqrt(dx*dx+dy*dy+dz*dz);
    dev = d - pdb_dist[i];
    dev2 = dev*dev;
    temp = -k_bnd*dev/d/(1.0-dev2/R0sq);

    fx = temp*dx;
    fy = temp*dy;
    fz = temp*dz;    

    force[ibead].x -= fx;
    force[ibead].y -= fy;
    force[ibead].z -= fz;

    force[jbead].x += fx;
    force[jbead].y += fy;
    force[jbead].z += fz;

  }
  
}

void random_force() {

  using namespace std;
  
  double var;
  int problem;

  var = sqrt(2.0*T*zeta/h);
  
  for( int i=1; i<=nbead; i++ ) {
    force[i].x += var*generator.gasdev();
    force[i].y += var*generator.gasdev();
    force[i].z += var*generator.gasdev();
    
  }
  
}

double rnd(double x)
{
  
  using namespace std;

  return ( (x>0) ? floor(x+0.5) : ceil(x-0.5) );

}

void update_neighbor_list() {

  double dx, dy, dz;
  double d2;
  int ibead, jbead, itype, jtype;
  double rcut, rcut2;

  nnl_att = 0;
  nnl_rep = 0;

  for (int i=1; i<=ncon_att; i++) {

    ibead = ibead_lj_nat[i];
    jbead = jbead_lj_nat[i];
    itype = itype_lj_nat[i];
    jtype = jtype_lj_nat[i];

    dx = unc_pos[jbead].x - unc_pos[ibead].x;
    dy = unc_pos[jbead].y - unc_pos[ibead].y;
    dz = unc_pos[jbead].z - unc_pos[ibead].z;

    dx -= boxl*rnd(dx/boxl);
    dy -= boxl*rnd(dy/boxl);
    dz -= boxl*rnd(dz/boxl);

    d2 = dx*dx+dy*dy+dz*dz;

    rcut = 3.2*lj_nat_pdb_dist[i];
    rcut2 = rcut*rcut;

    if (d2 < rcut2) {
      // add to neighbor list
      nnl_att++;
      ibead_neighbor_list_att[nnl_att] = ibead;
      jbead_neighbor_list_att[nnl_att] = jbead;
      itype_neighbor_list_att[nnl_att] = itype;
      jtype_neighbor_list_att[nnl_att] = jtype;
      nl_lj_nat_pdb_dist[nnl_att] = lj_nat_pdb_dist[i];
      nl_lj_nat_pdb_dist2[nnl_att] = lj_nat_pdb_dist2[i];
      nl_lj_nat_pdb_dist6[nnl_att] = lj_nat_pdb_dist6[i];
      nl_lj_nat_pdb_dist12[nnl_att] = lj_nat_pdb_dist12[i];
    }
  }

  for (int i=1; i<=ncon_rep; i++) {

    ibead = ibead_lj_non_nat[i];
    jbead = jbead_lj_non_nat[i];
    itype = itype_lj_non_nat[i];
    jtype = jtype_lj_non_nat[i];

    dx = unc_pos[jbead].x - unc_pos[ibead].x;
    dy = unc_pos[jbead].y - unc_pos[ibead].y;
    dz = unc_pos[jbead].z - unc_pos[ibead].z;

    dx -= boxl*rnd(dx/boxl);
    dy -= boxl*rnd(dy/boxl);
    dz -= boxl*rnd(dz/boxl);

    d2 = dx*dx+dy*dy+dz*dz;

    rcut = 3.2*sigma_rep[itype][jtype];
    rcut2 = rcut*rcut;

    if (d2 < rcut2) {
      // add to neighbor list
      nnl_rep++;
      ibead_neighbor_list_rep[nnl_rep] = ibead;
      jbead_neighbor_list_rep[nnl_rep] = jbead;
      itype_neighbor_list_rep[nnl_rep] = itype;
      jtype_neighbor_list_rep[nnl_rep] = jtype;
    }

  }
}

void update_cell_list() {

  double imcx, imcy, imcz, jmcx, jmcy, jmcz;
  double cdistx, cdisty, cdistz;

  int ibead, jbead, itype, jtype;

  nnl_att = 0;
  nnl_rep = 0;

  for (int i=1; i<=ncon_att; i++) {

    ibead = ibead_lj_nat[i];
    jbead = jbead_lj_nat[i];
    itype = itype_lj_nat[i];
    jtype = jtype_lj_nat[i];

    // get vector index from coordinates
    imcx = floor(unc_pos[ibead].x / lcell);
    imcy = floor(unc_pos[ibead].y / lcell);
    imcz = floor(unc_pos[ibead].z / lcell);

    // get vector index from coordinates
    jmcx = floor(unc_pos[jbead].x / lcell);
    jmcy = floor(unc_pos[jbead].y / lcell);
    jmcz = floor(unc_pos[jbead].z / lcell);

    cdistx = imcx - jmcx;
    cdisty = imcy - jmcy;
    cdistz = imcz - jmcz;

    cdistx -= ncell*rnd(cdistx/ncell);
    cdisty -= ncell*rnd(cdisty/ncell);
    cdistz -= ncell*rnd(cdistz/ncell);

    // check whether beads are in neighboring cells
    if (cdistx >= -1.0 && cdistx <= 1.0 && 
	cdisty >= -1.0 && cdisty <= 1.0 &&
	cdistz >= -1.0 && cdistz <= 1.0) {

      //  add ibead and jbead to neighbor list
      nnl_att++;
      ibead_neighbor_list_att[nnl_att] = ibead;
      jbead_neighbor_list_att[nnl_att] = jbead;
      itype_neighbor_list_att[nnl_att] = itype;
      jtype_neighbor_list_att[nnl_att] = jtype;
      nl_lj_nat_pdb_dist[nnl_att] = lj_nat_pdb_dist[i];
      nl_lj_nat_pdb_dist2[nnl_att] = lj_nat_pdb_dist2[i];
      nl_lj_nat_pdb_dist6[nnl_att] = lj_nat_pdb_dist6[i];
      nl_lj_nat_pdb_dist12[nnl_att] = lj_nat_pdb_dist12[i];
    }
  }

  for (int i=1; i<=ncon_rep; i++) {

    ibead = ibead_lj_non_nat[i];
    jbead = jbead_lj_non_nat[i];
    itype = itype_lj_non_nat[i];
    jtype = jtype_lj_non_nat[i];

    // get vector index from coordinates
    imcx = floor(unc_pos[ibead].x / lcell);
    imcy = floor(unc_pos[ibead].y / lcell);
    imcz = floor(unc_pos[ibead].z / lcell);

    // get vector index from coordinates
    jmcx = floor(unc_pos[jbead].x / lcell);
    jmcy = floor(unc_pos[jbead].y / lcell);
    jmcz = floor(unc_pos[jbead].z / lcell);

    cdistx = imcx - jmcx;
    cdisty = imcy - jmcy;
    cdistz = imcz - jmcz;

    // apply minimum image convention so cells wrap
    cdistx -= ncell*rnd(cdistx/ncell);
    cdisty -= ncell*rnd(cdisty/ncell);
    cdistz -= ncell*rnd(cdistz/ncell);

    // check whether beads are in neighboring cells
    if (cdistx >= -1.0 && cdistx <= 1.0 && 
	cdisty >= -1.0 && cdisty <= 1.0 &&
	cdistz >= -1.0 && cdistz <= 1.0) {
	  
      //  add ibead and jbead to neighbor list
      nnl_rep++;
      ibead_neighbor_list_rep[nnl_rep] = ibead;
      jbead_neighbor_list_rep[nnl_rep] = jbead;
      itype_neighbor_list_rep[nnl_rep] = itype;
      jtype_neighbor_list_rep[nnl_rep] = jtype;
    }
  }
}

void update_pair_list() {

  using namespace std;

  // declare host variables
  double dx, dy, dz;
  double d2;
  unsigned int ibead, jbead, itype, jtype;
  double rcut, rcut2;

  nil_att = 0;
  nil_rep = 0;

  // declare device variables

  // should be native distance
  for (int i=1; i<=nnl_att; i++) {

    ibead = ibead_neighbor_list_att[i];
    jbead = jbead_neighbor_list_att[i];
    itype = itype_neighbor_list_att[i];
    jtype = jtype_neighbor_list_att[i];

    dx = unc_pos[jbead].x - unc_pos[ibead].x;
    dy = unc_pos[jbead].y - unc_pos[ibead].y;
    dz = unc_pos[jbead].z - unc_pos[ibead].z;

    dx -= boxl*rnd(dx/boxl);
    dy -= boxl*rnd(dy/boxl);
    dz -= boxl*rnd(dz/boxl);

    d2 = dx*dx+dy*dy+dz*dz;

    rcut = 2.5*nl_lj_nat_pdb_dist[i];
    rcut2 = rcut*rcut;

    if (d2 < rcut2) {
      // add to interaction pair list
      nil_att++;
      ibead_pair_list_att[nil_att] = ibead;
      jbead_pair_list_att[nil_att] = jbead;
      itype_pair_list_att[nil_att] = itype;
      jtype_pair_list_att[nil_att] = jtype;
      pl_lj_nat_pdb_dist[nil_att] = nl_lj_nat_pdb_dist[i];
      pl_lj_nat_pdb_dist2[nil_att] = nl_lj_nat_pdb_dist2[i];
      pl_lj_nat_pdb_dist6[nil_att] = nl_lj_nat_pdb_dist6[i];
      pl_lj_nat_pdb_dist12[nil_att] = nl_lj_nat_pdb_dist12[i];
    }

  }

  for (int i=1; i<=nnl_rep; i++) {

    ibead = ibead_neighbor_list_rep[i];
    jbead = jbead_neighbor_list_rep[i];
    itype = itype_neighbor_list_rep[i];
    jtype = jtype_neighbor_list_rep[i];

    dx = unc_pos[jbead].x - unc_pos[ibead].x;
    dy = unc_pos[jbead].y - unc_pos[ibead].y;
    dz = unc_pos[jbead].z - unc_pos[ibead].z;

    dx -= boxl*rnd(dx/boxl);
    dy -= boxl*rnd(dy/boxl);
    dz -= boxl*rnd(dz/boxl);

    d2 = dx*dx+dy*dy+dz*dz;

    rcut = 2.5*sigma_rep[itype][jtype];
    rcut2 = rcut*rcut;

    if (d2 < rcut2) {
      // add to interaction pair list
      nil_rep++;
      ibead_pair_list_rep[nil_rep] = ibead;
      jbead_pair_list_rep[nil_rep] = jbead;
      itype_pair_list_rep[nil_rep] = itype;
      jtype_pair_list_rep[nil_rep] = jtype;
    }
  }
}

