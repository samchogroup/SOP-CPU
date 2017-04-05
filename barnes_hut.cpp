#include <iostream>
#include <cstring>
#include <math.h>
#include <cstdlib>
#include <cstdio>
#include <unistd.h>
#include "barnes_hut.h"

coord * find_subtree_coord(double width, coord *center, int bead_index){
  double half = width/2.0;
  coord * subtree = new coord();

  // Does not handle coordinates outside the box space
  subtree->x = unc_pos[bead_index].x > center->x - half && unc_pos[bead_index].x < center->x ? 0 : 1;
  subtree->y = unc_pos[bead_index].y > center->y - half && unc_pos[bead_index].y < center->y ? 0 : 1;
  subtree->z = unc_pos[bead_index].z > center->z - half && unc_pos[bead_index].z < center->z ? 0 : 1;

  // std::cout << "For bead " << bead_index << " coordinates are: " ;
  // std::cout << subtree->x << " " << subtree->y << " " << subtree->z << '\n';

  return subtree;
}

int find_subtree_index(coord *center, int current_node){

  int tree_child = (center->x * 4) + (center->y * 2) + center->z;
  tree_child++;

  int subtree_index = 8*current_node + tree_child;

  return subtree_index;
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

/*
  tree_index = index of the current subtree in the flattened array structures
  bead_index = index of the bead in the list of bead data
  box_center = tridimensional coordinate that marks the center of the box associated to the subtree
  width = widht of the box associated with the subtree
*/
void insert_bead_bhtree(int tree_index, int bead_index, coord *boxCenter, double width){
  // std::cout << "INSERTING BEAD " << bead_index << " ";
  // unc_pos[bead_index].print(); std::cout << '\n' << "In block: ";
  // boxCenter->print(); std::cout << " " << width << '\n';

  if (indices_bhtree[tree_index] == empty_node){
    // insert node here
    std::cout << "found an empty spot in " << tree_index << "... inserting" << '\n';
    indices_bhtree[tree_index] = bead_index;
    octet_count_bhtree[tree_index] = 1;
    octet_center_mass[tree_index].x = boxCenter->x;
    octet_center_mass[tree_index].y = boxCenter->y;
    octet_center_mass[tree_index].z = boxCenter->z;

  } else if (indices_bhtree[tree_index] == internal_node){
    std::cout << "found an internal node in " << tree_index << "... inserting recursively" << '\n';
    octet_count_bhtree[tree_index] = octet_count_bhtree[tree_index]+1;
    update_center_mass(tree_index, bead_index);

    coord *subtree = find_subtree_coord(width, boxCenter, bead_index);
    int subtree_index = find_subtree_index(subtree, tree_index);
    coord *newCenter = get_udpated_center(boxCenter, subtree, width);
    insert_bead_bhtree(subtree_index, bead_index, newCenter, width/2.0);
    delete newCenter;

  } else {
    std::cout << "found a leaf in " << tree_index << "... will make internal and insert recursively" << '\n';
    int a_bead_index = indices_bhtree[tree_index];
    indices_bhtree[tree_index] = internal_node;
    octet_count_bhtree[tree_index] = 2;
    update_center_mass(tree_index, bead_index);

    coord *subtreeA = find_subtree_coord(width, boxCenter, a_bead_index);
    coord *subtreeB = find_subtree_coord(width, boxCenter, bead_index);
    int subtree_indexA = find_subtree_index(subtreeA, tree_index);
    int subtree_indexB = find_subtree_index(subtreeB, tree_index);
    coord *newCenterA = get_udpated_center(boxCenter, subtreeA, width);
    coord *newCenterB = get_udpated_center(boxCenter, subtreeB, width);
    insert_bead_bhtree(subtree_indexA, a_bead_index, newCenterA, width/2.0);
    insert_bead_bhtree(subtree_indexB, bead_index, newCenterB, width/2.0);
    delete newCenterA; delete newCenterB;
  }
}

double get_initial_width(){
  double width = 0.0;
  float x,y,z;
  for (int i = 0; i < nbead; i++) {
    x = abs(unc_pos[i].x);
    y = abs(unc_pos[i].y);
    z = abs(unc_pos[i].z);
    width = x > width ? x : width;
    width = y > width ? y : width;
    width = z > width ? z : width;
  }
  return ceil(2*width);
}

void build_bh_tree(){

  int tree_index = 0; // all insertions start in the root of the tree
  double boxWidth = get_initial_width();
  coord *boxCenter = new coord(); // initial box center coordinate - 0, 0, 0
  boxCenter->x = 0.0;
  boxCenter->y = 0.0;
  boxCenter->z = 0.0;

  for (int i = 1; i < nbead; i++) {
    insert_bead_bhtree(tree_index, i, boxCenter, boxWidth);
    std::cin.get();
  }

}
