#ifndef BH_H
#define BH_H

#include "global.h"

extern int next_tree_index;
void build_bh_tree();
void insert_bead_bhtree(int tree_index, int bead_index, coord *boxCenter, int width);
void bh_update_pair_list();

#endif
