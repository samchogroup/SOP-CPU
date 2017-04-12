#include <cstdlib>
#include <math.h>
#include <iostream>
#include <cstdio>
#include "global.h"
#include "energy.h"

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

void clear_forces()
{

  using namespace std;

  for( int i=1; i<=nbead; i++ ) {
    force[i].x = 0.0;
    force[i].y = 0.0;
    force[i].z = 0.0;
  }

}

void set_potential()
{

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

void vdw_bh_energy()
{

  using namespace std;

  int ibead,jbead;
  int itype,jtype;
  double dx,dy,dz,d,d2,d6,d12;

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

    d2 = att_pl_bh_d2[i];
    d6 = att_pl_bh_d6[i];
    d12 = att_pl_bh_d12[i];

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

    d2 = rep_pl_bh_d2[i];
    d6 = rep_pl_bh_d6[i];
    d12 = rep_pl_bh_d12[i];

    e_vdw_rr_rep += coeff_rep[itype][jtype] * (sigma_rep12[itype][jtype]/d12+sigma_rep6[itype][jtype]/d6);

  }

  e_vdw_rr = e_vdw_rr_att + e_vdw_rr_rep;

  return;

}

void vdw_energy()
{
  using namespace std;

  // if(barnesHut){
  //   vdw_bh_energy();
  //   return;
  // }

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

void vdw_bh_forces()
{
  using namespace std;

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

    d2 = att_pl_bh_d2[i];
    d6 = att_pl_bh_d6[i];
    d12 = att_pl_bh_d12[i];

    rep_tol = sigma_rep2[itype][jtype]*tol;

    if( d2 < tol*lj_nat_pdb_dist2[i] ) continue;

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

    d2 = rep_pl_bh_d2[i];
    d6 = rep_pl_bh_d6[i];
    d12 = rep_pl_bh_d12[i];

    if( d2 <  rep_tol ) continue;

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

void vdw_forces()
{
  using namespace std;

  // if(barnesHut){
  //   vdw_bh_forces();
  //   return;
  // }

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
