#include <iostream>
#include <cstring>
#include <math.h>
#include <cstdlib>
#include <cstdio>
#include <unistd.h>
#include "barnes_hut.h"

int find_subtree_index(int width, coord *center, int bead_index, int current_node){
  int half = width/2;
  int xval, yval, zval;

  if (unc_pos[bead_index].x > center->x - half && unc_pos[bead_index].x < center->x) {
    xval = 0;
  } else if (unc_pos[bead_index].x > center->x && unc_pos[bead_index].x < center->x + half){
    xval = 1;
  } else {
    std::cout << "ERROR: OUT OF BOX" << '\n';
  }

  if (unc_pos[bead_index].y > center->y - half && unc_pos[bead_index].y < center->y) {
    yval = 0;
  } else if (unc_pos[bead_index].y > center->y && unc_pos[bead_index].y < center->y + half){
    yval = 1;
  } else {
    std::cout << "ERROR: OUT OF BOX" << '\n';
  }

  if (unc_pos[bead_index].z > center->z - half && unc_pos[bead_index].z < center->z) {
    zval = 0;
  } else if (unc_pos[bead_index].z > center->z && unc_pos[bead_index].z < center->z + half){
    zval = 1;
  } else {
    std::cout << "ERROR: OUT OF BOX" << '\n';
  }

  int tree_child = (xval * 4) + (yval * 2) + zval;
  tree_child++;

  int subtree_index = 8*current_node + tree_child;

  return subtree_index;

}

void insert_bead_bhtree(int tree_index, int bead_index, coord *boxCenter, int width){
  if (indices_bhtree[tree_index] == empty_node){
    // insert node here
    std::cout << "found an empty spot in " << tree_index << "... inserting" << '\n';
    indices_bhtree[tree_index] = bead_index;
    octet_count_bhtree[tree_index] = 1;
    octet_width_bhtree[tree_index] = width;
    octet_avg_pos[tree_index].x = boxCenter->x;
    octet_avg_pos[tree_index].y = boxCenter->y;
    octet_avg_pos[tree_index].z = boxCenter->z;

  } else if (indices_bhtree[tree_index] == inner_node){
    std::cout << "found an inner node in " << tree_index << "... inserting recursively" << '\n';
    octet_count_bhtree[tree_index] = octet_count_bhtree[tree_index]+1;
    // update center of mass
    // find subtree
    // insert


  } else {
    std::cout << "found a leaf in " << tree_index << "... will make inner and insert recursively" << '\n';
    indices_bhtree[tree_index] = inner_node;
    octet_count_bhtree[tree_index] = 2;
    // update center of mass
    // find subtree for both
    // insert both
  }
}

// indices_bhtree = new int[2*nbead];
// octet_count_bhtree = new double[2*nbead];
// octet_width_bhtree = new double[2*nbead];
// octet_avg_pos = new coord[2*nbead];

void build_bh_tree(){

  int tree_index = 0; // all insertions start in the root of the tree
  int boxWidth = 512; // the initial box width
  coord *boxCenter = new coord(); // initial box center coordinate - 0, 0, 0
  boxCenter->x = 0.0;
  boxCenter->y = 0.0;
  boxCenter->z = 0.0;

  for (int i = 0; i < nbead; i++) {
    insert_bead_bhtree(tree_index, i, boxCenter, boxWidth);
  }

}
