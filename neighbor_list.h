#ifndef NLIST_H
#define NLIST_H

#include "global.h"


void update_neighbor_list();
void update_neighbor_list_gpu();

void ksScanAuxExc (int *X, int *Y, int InputSize, int *S);
void ksScanAuxInc (int *X, int *Y, int InputSize, int *S);
void ksScanExc (int *X, int *Y, int InputSize);
void ksScanInc (int *X, int *Y, int InputSize);
void sumIt (int *Y, int *S, int InputSize);

void copyElements(int *dev_index, int *dev_value, int *dev_output, int *dev_result, int N);
void copyElements(double *dev_index, int *dev_value, int *dev_output, double *dev_result, int N);

void calculate_array_native(int *ibead_lj_nat, int *jbead_lj_nat, int *itype_lj_nat, int *jtype_lj_nat, coord *unc_pos, double *lj_nat_pdb_dist, int *value, int boxl, int N);
void array_native_kernel(int *dev_ibead_lj_nat, int *dev_jbead_lj_nat, int *dev_itype_lj_nat, int *dev_jtype_lj_nat, coord *dev_unc_pos, double *dev_lj_nat_pdb_dist, int *dev_value, int boxl, int N);

void calculate_array_non_native(int *ibead_lj_nat, int *jbead_lj_nat, int *itype_lj_nat, int *jtype_lj_nat, coord *unc_pos, int *value, int boxl, int N);
void array_non_native_kernel(int *dev_ibead_lj_nat, int *dev_jbead_lj_nat, int *dev_itype_lj_nat, int *dev_jtype_lj_nat, coord *dev_unc_pos, int *dev_value, int boxl, int N);

void hier_ks_scan(int *dev_X, int *dev_Y, int N, int re);
int compact(int *index, int *value, int N, int *&result);
int compact(double *index, int *value, int N, double *&result);

void dummy(int *dev_ibead_lj_nat, int *dev_jbead_lj_nat, int *dev_itype_lj_nat, int *dev_jtype_lj_nat, coord *dev_unc_pos, double *dev_lj_nat_pdb_dist, int *&dev_value, int boxl, int N, int nbead);

#endif
