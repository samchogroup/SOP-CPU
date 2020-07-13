#ifndef NLIST_H
#define NLIST_H

#include "global.h"


void update_neighbor_list();
void update_neighbor_list_gpu();

__global__ void ksScanAuxExc (int *X, int *Y, int InputSize, int *S);
__global__ void ksScanAuxInc (int *X, int *Y, int InputSize, int *S);
__global__ void ksScanExc (int *X, int *Y, int InputSize);
__global__ void ksScanInc (int *X, int *Y, int InputSize);
__global__ void sumIt (int *Y, int *S, int InputSize);

__global__ void copyElements(int *dev_index, int *dev_value, int *dev_output, int *dev_result, int N);
__global__ void copyElements(double *dev_index, int *dev_value, int *dev_output, double *dev_result, int N);

void calculate_array_native(int *ibead_lj_nat, int *jbead_lj_nat, int *itype_lj_nat, int *jtype_lj_nat, float3 *unc_pos, double *lj_nat_pdb_dist, int *value, int boxl, int N);
__global__ void array_native_kernel(int *dev_ibead_lj_nat, int *dev_jbead_lj_nat, int *dev_itype_lj_nat, int *dev_jtype_lj_nat, float3 *dev_unc_pos, double *dev_lj_nat_pdb_dist, int *dev_value, int boxl, int N);

void calculate_array_non_native(int *ibead_lj_nat, int *jbead_lj_nat, int *itype_lj_nat, int *jtype_lj_nat, float3 *unc_pos, int *value, int boxl, int N);
__global__ void array_non_native_kernel(int *dev_ibead_lj_nat, int *dev_jbead_lj_nat, int *dev_itype_lj_nat, int *dev_jtype_lj_nat, float3 *dev_unc_pos, int *dev_value, int boxl, int N);

void hier_ks_scan(int *dev_X, int *dev_Y, int N, int re);
int compact(int *index, int *value, int N, int *&result);
int compact(double *index, int *value, int N, double *&result);

__global__ void dummy(int *dev_ibead_lj_nat, int *dev_jbead_lj_nat, int *dev_itype_lj_nat, int *dev_jtype_lj_nat, float3 *dev_unc_pos, double *dev_lj_nat_pdb_dist, int *&dev_value, int boxl, int N, int nbead);

#endif
