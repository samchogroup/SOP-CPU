#include <iostream>
#include <cstring>
#include <math.h>
#include <cstdlib>
#include <cstdio>
#include <unistd.h>
#include "barnes_hut.h"

int next_tree_index = 0;
double rootWidth = 0;

int get_next_tree_index(){
  int old = next_tree_index + 1;
  next_tree_index += 8; //creates space for 8 new cells
  return old;
}

coord * find_subtree_coord(double width, coord *center, int bead_index){
  coord * subtree = new coord();

  subtree->x = unc_pos[bead_index].x < center->x ? 0 : 1;
  subtree->y = unc_pos[bead_index].y < center->y ? 0 : 1;
  subtree->z = unc_pos[bead_index].z < center->z ? 0 : 1;

  return subtree;
}

int get_subtree_index(coord *sub_block, int base_index){
  int tree_child = base_index + (sub_block->x * 4) + (sub_block->y * 2) + sub_block->z;
  return tree_child;
}

void update_center_mass(int current_node, int bead_index){
  octet_center_mass[current_node].x += unc_pos[bead_index].x;
  octet_center_mass[current_node].y += unc_pos[bead_index].y;
  octet_center_mass[current_node].z += unc_pos[bead_index].z;
}

coord * get_udpated_center(coord *boxCenter, coord* subtree, double width){
  double quarter = width/4.0;
  coord * newCenter = new coord();

  newCenter->x = subtree->x == 1.0 ? boxCenter->x + quarter : boxCenter->x - quarter;
  newCenter->y = subtree->y == 1.0 ? boxCenter->y + quarter : boxCenter->y - quarter;
  newCenter->z = subtree->z == 1.0 ? boxCenter->z + quarter : boxCenter->z - quarter;

  return newCenter;
}

void insert_bead_bhtree(int tree_index, int bead_index, coord *boxCenter, double width){
  if (indices_bhtree[tree_index] == empty_cell){
    /* Inseart bead with negative index */
    indices_bhtree[tree_index] = -bead_index;
    octet_count_bhtree[tree_index] = 1;
    octet_center_mass[tree_index].x = boxCenter->x;
    octet_center_mass[tree_index].y = boxCenter->y;
    octet_center_mass[tree_index].z = boxCenter->z;

  } else if (indices_bhtree[tree_index] < 0) {
    /* Index is negative if there's a node in the cell, new positive index will be added for new cell */
    int a_bead_index = -indices_bhtree[tree_index];
    indices_bhtree[tree_index] = get_next_tree_index();
    octet_count_bhtree[tree_index] = 2;
    update_center_mass(tree_index, bead_index);

    coord *subtreeA = find_subtree_coord(width, boxCenter, a_bead_index);
    coord *subtreeB = find_subtree_coord(width, boxCenter, bead_index);
    int subtree_indexA = get_subtree_index(subtreeA, indices_bhtree[tree_index]);
    int subtree_indexB = get_subtree_index(subtreeB, indices_bhtree[tree_index]);
    coord *newCenterA = get_udpated_center(boxCenter, subtreeA, width);
    coord *newCenterB = get_udpated_center(boxCenter, subtreeB, width);
    insert_bead_bhtree(subtree_indexA, a_bead_index, newCenterA, width/2.0);
    insert_bead_bhtree(subtree_indexB, bead_index, newCenterB, width/2.0);
    delete newCenterA; delete newCenterB;

  } else {
    /* Positive index, internal node found */
    octet_count_bhtree[tree_index] = octet_count_bhtree[tree_index]+1;
    update_center_mass(tree_index, bead_index);

    coord *subtree = find_subtree_coord(width, boxCenter, bead_index);
    int subtree_index = get_subtree_index(subtree, indices_bhtree[tree_index]);
    coord *newCenter = get_udpated_center(boxCenter, subtree, width);
    insert_bead_bhtree(subtree_index, bead_index, newCenter, width/2.0);
    delete newCenter;
  }
}

double set_initial_width(){
  double width = 0.0;
  double x,y,z;
  for (int i = 0; i <= nbead; i++) {
    x = fabs(unc_pos[i].x);
    y = fabs(unc_pos[i].y);
    z = fabs(unc_pos[i].z);
    width = x > width ? x : width;
    width = y > width ? y : width;
    width = z > width ? z : width;
  }
  rootWidth = ceil(2*width);
}

void build_bh_tree(){
  reinserted++;
  /* reset tree */
  std::fill_n(indices_bhtree, 16*nbead, empty_cell);
  next_tree_index = 0;
  int tree_index = 0;
  set_initial_width();
  coord *boxCenter = new coord();
  boxCenter->x = 0.0;
  boxCenter->y = 0.0;
  boxCenter->z = 0.0;

  for (int i = 1; i <= nbead; i++) {
    insert_bead_bhtree(tree_index, i, boxCenter, rootWidth);
  }

  rebuild = 0;
  delete boxCenter;
}

coord * calculateDistances(double d2){
  coord * distances = new coord();

  double d6 = d2*d2*d2;
  double d12 = d6*d6;

  distances->x = d2;
  distances->y = d6;
  distances->z = d12;

  return distances;
}

double getDistanceSquared(int ibead, int jbead){
  double dx = unc_pos[jbead].x - unc_pos[ibead].x;
  double dy = unc_pos[jbead].y - unc_pos[ibead].y;
  double dz = unc_pos[jbead].z - unc_pos[ibead].z;

  dx -= boxl*rnd(dx/boxl);
  dy -= boxl*rnd(dy/boxl);
  dz -= boxl*rnd(dz/boxl);

  double d2 = dx*dx+dy*dy+dz*dz;
  return d2;
}

coord * getPreCalculatedDist(int tree_index){
  coord * distances = new coord();

  distances->x = aux_tree_d2[tree_index];
  distances->y = aux_tree_d6[tree_index];
  distances->z = aux_tree_d12[tree_index];

  return distances;
}

coord * getDistances(int tree_index, int ibead, int jbead, coord *boxCenter, double width){

  if (indices_bhtree[tree_index] == empty_cell){
    /* bead moved out of cell, mark build for rebuild and return distance */
    rebuild = 1;
    return calculateDistances(getDistanceSquared(ibead,jbead));
  }

 if (indices_bhtree[tree_index] < empty_cell) {
   /* leaf found, calculate distance and return */

    if (jbead != -indices_bhtree[tree_index]) {
      /* bead moved out of cell, mark build for rebuild and return distance*/
      rebuild = 1;
      return calculateDistances(getDistanceSquared(ibead,jbead));
    }

    /* found correct leaf */
    individual++;
    return calculateDistances(getDistanceSquared(ibead,jbead));

  } else {
    /* internal node (cell) found */

    if (aux_tree_d2[tree_index] != 0) {
      /* already calculated */
      approximated++;
      return getPreCalculatedDist(tree_index);
    }

    /* needs calculation */
    double s2 = width*width;

    double dx = (octet_center_mass[tree_index].x/octet_count_bhtree[tree_index]) - unc_pos[ibead].x;
    double dy = (octet_center_mass[tree_index].y/octet_count_bhtree[tree_index]) - unc_pos[ibead].y;
    double dz = (octet_center_mass[tree_index].z/octet_count_bhtree[tree_index]) - unc_pos[ibead].z;

    dx -= boxl*rnd(dx/boxl);
    dy -= boxl*rnd(dy/boxl);
    dz -= boxl*rnd(dz/boxl);

    double d2 = dx*dx+dy*dy+dz*dz;

    if(s2/d2 < theta2){
      /* cell is far enough to use this distance for calculations */
      double d6 = d2*d2*d2;
      double d12 = d6*d6;
      aux_tree_d2[tree_index] = d2;
      aux_tree_d6[tree_index] = d6;
      aux_tree_d12[tree_index] = d12;

      approximated++;
      return getPreCalculatedDist(tree_index);

    } else {
      /* cell is close, search deeper */
      coord *subtree = find_subtree_coord(width, boxCenter, jbead);
      int subtree_index = get_subtree_index(subtree, indices_bhtree[tree_index]);
      coord *newCenter = get_udpated_center(boxCenter, subtree, width);
      coord *distances = getDistances(subtree_index, ibead, jbead, newCenter, width/2);
      delete newCenter;
      return distances;
    }
  }
}

void bh_update_pair_list(){
  std::fill_n(aux_tree_d2, 16*nbead, empty_cell);
  nil_att = 0;
  nil_rep = 0;
  unsigned int ibead, jbead, itype, jtype;

  coord *boxCenter = new coord();
  boxCenter->x = 0.0;
  boxCenter->y = 0.0;
  boxCenter->z = 0.0;

  /* attractive interactions */

  for (int i=1; i<=ncon_att; i++) {
    if (ibead != ibead_lj_nat[i]) {
      /* the group of interactions ibead-x is done, aux tree with distances is zeroed out */
      std::fill_n(aux_tree_d2, 16*nbead, empty_cell);
    }

    ibead = ibead_lj_nat[i];
    jbead = jbead_lj_nat[i];
    itype = itype_lj_nat[i];
    jtype = jtype_lj_nat[i];

    coord *distances = getDistances(0, ibead, jbead, boxCenter, rootWidth);

    nil_att++;
    ibead_pair_list_att[nil_att] = ibead;
    jbead_pair_list_att[nil_att] = jbead;
    itype_pair_list_att[nil_att] = itype;
    jtype_pair_list_att[nil_att] = jtype;
    pl_lj_nat_pdb_dist[nil_att] = lj_nat_pdb_dist[i];
    pl_lj_nat_pdb_dist2[nil_att] = lj_nat_pdb_dist2[i];
    pl_lj_nat_pdb_dist6[nil_att] = lj_nat_pdb_dist6[i];
    pl_lj_nat_pdb_dist12[nil_att] = lj_nat_pdb_dist12[i];

    att_pl_bh_d2[nil_att] = distances->x;
    att_pl_bh_d6[nil_att] = distances->y;
    att_pl_bh_d12[nil_att] = distances->z;

    delete distances;
  }

  /* aux tree reset for repulsive pair list */
  std::fill_n(aux_tree_d2, 16*nbead, empty_cell);

  /* repulsive interactions */
  for (int i=1; i<=ncon_rep; i++) {

    if (ibead != ibead_lj_nat[i]) {
      /* the group of interactions ibead-x is done, aux tree with distances is zeroed out */
      std::fill_n(aux_tree_d2, 16*nbead, empty_cell);
    }

    ibead = ibead_lj_non_nat[i];
    jbead = jbead_lj_non_nat[i];
    itype = itype_lj_non_nat[i];
    jtype = jtype_lj_non_nat[i];

    coord *distances = getDistances(0, ibead, jbead, boxCenter, rootWidth);

    //rcut?
    nil_rep++;
    ibead_pair_list_rep[nil_rep] = ibead;
    jbead_pair_list_rep[nil_rep] = jbead;
    itype_pair_list_rep[nil_rep] = itype;
    jtype_pair_list_rep[nil_rep] = jtype;

    rep_pl_bh_d2[nil_rep] = distances->x;
    rep_pl_bh_d6[nil_rep] = distances->y;
    rep_pl_bh_d12[nil_rep] = distances->z;

    delete distances;

  }

  delete boxCenter;
}
