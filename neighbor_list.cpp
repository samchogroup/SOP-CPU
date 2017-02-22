#include <math.h>
#include <cstdlib>
#include "neighbor_list.h"
#include "global.h"

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
