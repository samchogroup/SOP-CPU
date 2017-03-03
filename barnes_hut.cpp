#include <cstring>
#include <cstdlib>
#include <cstdio>
#include "barnes_hut.h"

// indices_bhtree = new int[2*nbead];
// octet_count_bhtree = new double[2*nbead];
// octet_width_bhtree = new double[2*nbead];
// octet_avg_pos = new coord[2*nbead];

void insert_bead_bhtree(int tree_index, int bead_index, coord *boxCenter, int width){
  if (tree_index == empty_node){
    // insert node here
  } else if (tree_index == inner_node){
    // update count, pos and recursively insert
  } else {
    // leaf, index = -2, width = width/2, count = 2, average pos, insert both recursively
  }
}

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
