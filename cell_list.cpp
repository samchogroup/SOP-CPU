#include <math.h>
#include <cstdlib>
#include "cell_list.h"
#include "global.h"

void update_cell_list() {

  double imcx, imcy, imcz, jmcx, jmcy, jmcz;
  double cdistx, cdisty, cdistz;

  int ibead, jbead, itype, jtype;

  nnl_att = 0;
  nnl_rep = 0;

  for (int i=0; i<ncon_att; i++) {

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

  for (int i=0; i<ncon_rep; i++) {

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
