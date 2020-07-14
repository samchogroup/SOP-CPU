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

  // calculations for native (attractiction) contacts
  for (int i=0; i<ncon_att; i++) {
    // record sigma for ibead and jbead
    ibead = ibead_lj_nat[i];
    jbead = jbead_lj_nat[i];

    // record type of bead for ibead and jbead
    itype = itype_lj_nat[i];
    jtype = jtype_lj_nat[i];
    
    // calculate distance in x, y, and z for ibead and jbead
    dx = unc_pos[jbead-1].x - unc_pos[ibead-1].x;
    dy = unc_pos[jbead-1].y - unc_pos[ibead-1].y;
    dz = unc_pos[jbead-1].z - unc_pos[ibead-1].z;

    // apply periodic boundary conditions to dx, dy, and dz
    dx -= boxl*rnd(dx/boxl);
    dy -= boxl*rnd(dy/boxl);
    dz -= boxl*rnd(dz/boxl);

    // compute square of distance between ibead and jbead
    d2 = dx*dx+dy*dy+dz*dz;

    /* 
    Compute the cutoff distance for the given bead
    This is based off of lj_nat_pdb_dist[i], which is the distance 
    from ibead to jbead in the resulting folded structure
    */
    rcut = 3.2*lj_nat_pdb_dist[i];

    // square cutoff distance, since sqrt(d2) is computationally expensive
    rcut2 = rcut*rcut;

    // checks if distance squared is less than the cutoff distance squared
    if (d2 < rcut2) {
      // add pair to respective attraction neighbor lists
      ibead_neighbor_list_att[nnl_att] = ibead;
      jbead_neighbor_list_att[nnl_att] = jbead;
      
      // record type of each bead
      itype_neighbor_list_att[nnl_att] = itype;
      jtype_neighbor_list_att[nnl_att] = jtype;

      // record values, so that calculatons are not repeated (look-up table)
      nl_lj_nat_pdb_dist[nnl_att] = lj_nat_pdb_dist[i];
      nl_lj_nat_pdb_dist2[nnl_att] = lj_nat_pdb_dist2[i];
      nl_lj_nat_pdb_dist6[nnl_att] = lj_nat_pdb_dist6[i];
      nl_lj_nat_pdb_dist12[nnl_att] = lj_nat_pdb_dist12[i];
      // add to neighbor list
      nnl_att++;
    }
  }

  // calculations for non-native (repulsive) contacts
  for (int i=0; i<ncon_rep; i++) {
    // record sigma for ibead and jbead
    ibead = ibead_lj_non_nat[i];
    jbead = jbead_lj_non_nat[i];

    // record type of bead for ibead and jbead
    itype = itype_lj_non_nat[i];
    jtype = jtype_lj_non_nat[i];

    // calculate distance in x, y, and z for ibead and jbead
    dx = unc_pos[jbead-1].x - unc_pos[ibead-1].x;
    dy = unc_pos[jbead-1].y - unc_pos[ibead-1].y;
    dz = unc_pos[jbead-1].z - unc_pos[ibead-1].z;

    // apply periodic boundary conditions to dx, dy, and dz
    dx -= boxl*rnd(dx/boxl);
    dy -= boxl*rnd(dy/boxl);
    dz -= boxl*rnd(dz/boxl);

    // compute square of distance between ibead and jbead
    d2 = dx*dx+dy*dy+dz*dz;

    /* 
    Compute the cutoff distance for the given bead
    This is based off of sigma_rep[itype][jtype],
    is based on the sigma for the types of ibead and jbead
    */
    rcut = 3.2*sigma_rep[itype][jtype];

    // square cutoff distance, since sqrt(d2) is computationally expensive
    rcut2 = rcut*rcut;

    // checks if distance squared is less than the cutoff distance squared
    if (d2 < rcut2) {
      // add pair to respective repulsive neighbor lists
      ibead_neighbor_list_rep[nnl_rep] = ibead;
      jbead_neighbor_list_rep[nnl_rep] = jbead;

      // record type of each bead
      itype_neighbor_list_rep[nnl_rep] = itype;
      jtype_neighbor_list_rep[nnl_rep] = jtype;
      // add to neighbor list
      nnl_rep++;

    }

  }
}
