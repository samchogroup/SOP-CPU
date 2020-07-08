#include <math.h>
#include <cstdlib>
#include "neighbor_list.h"
#include "global.h"
#include <stdio.h>
#include <stdlib.h>

#define SECTION_SIZE 1024

void update_neighbor_list() {

  double dx, dy, dz;
  double d2;
  int ibead, jbead, itype, jtype;
  double rcut, rcut2;

  nnl_att = 0;
  nnl_rep = 0;

  // calculations for native (attractiction) contacts
  for (int i=1; i<=ncon_att; i++) {
    // record sigma for ibead and jbead
    ibead = ibead_lj_nat[i];
    jbead = jbead_lj_nat[i];

    // record type of bead for ibead and jbead
    itype = itype_lj_nat[i];
    jtype = jtype_lj_nat[i];
    
    // calculate distance in x, y, and z for ibead and jbead
    dx = unc_pos[jbead].x - unc_pos[ibead].x;
    dy = unc_pos[jbead].y - unc_pos[ibead].y;
    dz = unc_pos[jbead].z - unc_pos[ibead].z;

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
      // add to neighbor list
      nnl_att++;
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
    }
  }

  // calculations for non-native (repulsive) contacts
  for (int i=1; i<=ncon_rep; i++) {
    // record sigma for ibead and jbead
    ibead = ibead_lj_non_nat[i];
    jbead = jbead_lj_non_nat[i];

    // record type of bead for ibead and jbead
    itype = itype_lj_non_nat[i];
    jtype = jtype_lj_non_nat[i];

    // calculate distance in x, y, and z for ibead and jbead
    dx = unc_pos[jbead].x - unc_pos[ibead].x;
    dy = unc_pos[jbead].y - unc_pos[ibead].y;
    dz = unc_pos[jbead].z - unc_pos[ibead].z;

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
      // add to neighbor list
      nnl_rep++;

      // add pair to respective repulsive neighbor lists
      ibead_neighbor_list_rep[nnl_rep] = ibead;
      jbead_neighbor_list_rep[nnl_rep] = jbead;

      // record type of each bead
      itype_neighbor_list_rep[nnl_rep] = itype;
      jtype_neighbor_list_rep[nnl_rep] = jtype;
    }

  }
}

void update_neighbor_list_gpu() {
  int N;

  // Declare pointers for dev_output and dev_value arrays
  int *dev_value;
  int *dev_output;
  int *dev_unc_pos;

  // Start NL Update for Attractive Pairs
  int *dev_ibead_lj_nat;
  int *dev_jbead_lj_nat;
  int *dev_itype_lj_nat;
  int *dev_jtype_lj_nat;
  int *dev_lj_nat_pdb_dist;

  // Declare local variables

  // Calculate array size
  N = ncon_att+1;

  int size = (N) * sizeof(int);

  // Allocate dev_ arrays
  // TODO: Calculate correct size for each array
  cudaMalloc((void**)&dev_value, size);
  cudaMalloc((void**)&dev_output, size);
  cudaMalloc((void**)&dev_ibead_lj_nat, size);
  cudaMalloc((void**)&dev_jbead_lj_nat, size);
  cudaMalloc((void**)&dev_itype_lj_nat, size);
  cudaMalloc((void**)&dev_jtype_lj_nat, size);
  cudaMalloc((void**)&dev_unc_pos, size);
  cudaMalloc((void**)&dev_lj_nat_pdb_dist, size);

  // Copy arrays to dev_ arrays
  cudaMemcpy(dev_ibead_lj_nat, ibead_lj_nat, size, cudaMemcpyHostToDevice);
  cudaMemcpy(dev_jbead_lj_nat, jbead_lj_nat, size, cudaMemcpyHostToDevice);
  cudaMemcpy(dev_itype_lj_nat, itype_lj_nat, size, cudaMemcpyHostToDevice);
  cudaMemcpy(dev_jtype_lj_nat, jtype_lj_nat, size, cudaMemcpyHostToDevice);
  cudaMemcpy(dev_unc_pos, unc_pos, size, cudaMemcpyHostToDevice);
  cudaMemcpy(dev_lj_nat_pdb_dist, lj_nat_pdb_dist, size, cudaMemcpyHostToDevice);

  int threads = (int)min(N, SECTION_SIZE);
  int blocks = (int)ceil(1.0*N/SECTION_SIZE);
  calculate_array_native<<<threads, blocks>>>(dev_ibead_lj_nat, dev_jbead_lj_nat, dev_itype_lj_nat, dev_jtype_lj_nat, dev_unc_pos, dev_lj_nat_pdb_dist, 
                                              dev_value, boxl, N);

  // Free memory used for calculating binary values
  cudaFree(dev_unc_pos);
  

  // Run Scan on array
  hier_ks_scan(dev_value, dev_output, N, 0);
  int arrSize_att;
  int endVal_att;
  cudaMemcpy(&arrSize_att, &dev_output[N-1], sizeof(int), cudaMemcpyDeviceToHost); 
  cudaMemcpy(&endVal_att, &dev_value[N-1], sizeof(int), cudaMemcpyDeviceToHost); 

  // Increment arrSize by 1 if needed
  if(endVal_att){
      arrSize_att++;
  }

  // Declare more local variables
  int *dev_ibead_neighbor_list_att;
  int *dev_jbead_neighbor_list_att;
  int *dev_itype_neighbor_list_att;
  int *dev_jtype_neighbor_list_att;

  int *dev_lj_nat_pdb_dist2;
  int *dev_lj_nat_pdb_dist6;
  int *dev_lj_nat_pdb_dist12;

  // Allocate arrays
  cudaMalloc((void**)&dev_ibead_neighbor_list_att, size);
  cudaMalloc((void**)&dev_jbead_neighbor_list_att, size);
  cudaMalloc((void**)&dev_itype_neighbor_list_att, size);
  cudaMalloc((void**)&dev_jtype_neighbor_list_att, size);

  cudaMalloc((void**)&dev_lj_nat_pdb_dist2, size);
  cudaMalloc((void**)&dev_lj_nat_pdb_dist6, size);
  cudaMalloc((void**)&dev_lj_nat_pdb_dist12, size);

  // Copy to from host to device
  cudaMemcpy(dev_ibead_neighbor_list_att, ibead_neighbor_list_att, size, cudaMemcpyHostToDevice);
  cudaMemcpy(dev_jbead_neighbor_list_att, jbead_neighbor_list_att, size, cudaMemcpyHostToDevice);
  cudaMemcpy(dev_itype_neighbor_list_att, itype_neighbor_list_att, size, cudaMemcpyHostToDevice);
  cudaMemcpy(dev_jtype_neighbor_list_att, jtype_neighbor_list_att, size, cudaMemcpyHostToDevice);

  cudaMemcpy(dev_lj_nat_pdb_dist2, lj_nat_pdb_dist2, size, cudaMemcpyHostToDevice);
  cudaMemcpy(dev_lj_nat_pdb_dist6, lj_nat_pdb_dist6, size, cudaMemcpyHostToDevice);
  cudaMemcpy(dev_lj_nat_pdb_dist12, lj_nat_pdb_dist12, size, cudaMemcpyHostToDevice);

  // Copy elements
  threads = (int)min(N, SECTION_SIZE);
  blocks = (int)ceil(1.0*N/SECTION_SIZE);

  // Kernel to copy elements from dev_index to dev_output if their corresponding dev_value is 1
  copyElements<<<blocks, threads>>>(dev_ibead_lj_nat, dev_value, dev_output, dev_ibead_neighbor_list_att, N);
  copyElements<<<blocks, threads>>>(dev_jbead_lj_nat, dev_value, dev_output, dev_jbead_neighbor_list_att, N);
  copyElements<<<blocks, threads>>>(dev_jbead_lj_nat, dev_value, dev_output, dev_itype_neighbor_list_att, N);
  copyElements<<<blocks, threads>>>(dev_ibead_lj_nat, dev_value, dev_output, dev_jtype_neighbor_list_att, N);

  copyElements<<<blocks, threads>>>(dev_nl_lj_nat_pdb_dist, dev_value, dev_output, dev_lj_nat_pdb_dist, N);
  copyElements<<<blocks, threads>>>(dev_nl_lj_nat_pdb_dist2, dev_value, dev_output, dev_lj_nat_pdb_dist2, N);
  copyElements<<<blocks, threads>>>(dev_nl_lj_nat_pdb_dist6, dev_value, dev_output, dev_lj_nat_pdb_dist6, N);
  copyElements<<<blocks, threads>>>(dev_nl_lj_nat_pdb_dist12, dev_value, dev_output, dev_lj_nat_pdb_dist12, N);

  // Copy from device to host
  cudaMemcpy(ibead_neighbor_list_att, dev_ibead_neighbor_list_att, arrSize_att*sizeof(int), cudaMemcpyDeviceToHost);
  cudaMemcpy(jbead_neighbor_list_att, dev_jbead_neighbor_list_att, arrSize_att*sizeof(int), cudaMemcpyDeviceToHost);
  cudaMemcpy(itype_neighbor_list_att, dev_itype_neighbor_list_att, arrSize_att*sizeof(int), cudaMemcpyDeviceToHost);
  cudaMemcpy(jtype_neighbor_list_att, dev_jtype_neighbor_list_att, arrSize_att*sizeof(int), cudaMemcpyDeviceToHost);
  
  cudaMemcpy(lj_nat_pdb_dist, dev_lj_nat_pdb_dist, arrSize_att*sizeof(int), cudaMemcpyDeviceToHost);
  cudaMemcpy(lj_nat_pdb_dist2, dev_lj_nat_pdb_dist2, arrSize_att*sizeof(int), cudaMemcpyDeviceToHost);
  cudaMemcpy(lj_nat_pdb_dist6, dev_lj_nat_pdb_dist6, arrSize_att*sizeof(int), cudaMemcpyDeviceToHost);
  cudaMemcpy(lj_nat_pdb_dist12, dev_lj_nat_pdb_dist12, arrSize_att*sizeof(int), cudaMemcpyDeviceToHost);

  // Free the rest of the memory
  cudaFree(dev_ibead_lj_nat);
  cudaFree(dev_jbead_lj_nat);
  cudaFree(dev_itype_lj_nat);
  cudaFree(dev_jtype_lj_nat);
  cudaFree(dev_lj_nat_pdb_dist);
  cudaFree(dev_output);
  cudaFree(dev_value);

  /*************************************************************************************************************************************/

  // Start NL Update for Repulsive Pairs
  int *dev_ibead_lj_non_nat;
  int *dev_jbead_lj_non_nat;
  int *dev_itype_lj_non_nat;
  int *dev_jtype_lj_non_nat;
  int *dev_lj_non_nat_pdb_dist;

  // Calculate array size
  N = ncon_rep+1;

  int size = (N) * sizeof(int);

  // Allocate dev_ arrays
  // TODO: Calculate correct size for each array
  cudaMalloc((void**)&dev_value, size);
  cudaMalloc((void**)&dev_output, size);
  cudaMalloc((void**)&dev_ibead_lj_non_nat, size);
  cudaMalloc((void**)&dev_jbead_lj_non_nat, size);
  cudaMalloc((void**)&dev_itype_lj_non_nat, size);
  cudaMalloc((void**)&dev_jtype_lj_non_nat, size);
  cudaMalloc((void**)&dev_unc_pos, size);
  cudaMalloc((void**)&dev_lj_nat_pdb_dist, size);

  // Copy arrays to dev_ arrays
  cudaMemcpy(dev_ibead_lj_non_nat, ibead_lj_nat, size, cudaMemcpyHostToDevice);
  cudaMemcpy(dev_jbead_lj_non_nat, jbead_lj_nat, size, cudaMemcpyHostToDevice);
  cudaMemcpy(dev_itype_lj_non_nat, itype_lj_nat, size, cudaMemcpyHostToDevice);
  cudaMemcpy(dev_jtype_lj_non_nat, jtype_lj_nat, size, cudaMemcpyHostToDevice);
  cudaMemcpy(dev_unc_pos, unc_pos, size, cudaMemcpyHostToDevice);
  cudaMemcpy(dev_lj_non_nat_pdb_dist, lj_nat_pdb_dist, size, cudaMemcpyHostToDevice);

  threads = (int)min(N, SECTION_SIZE);
  blocks = (int)ceil(1.0*N/SECTION_SIZE);
  calculate_array_non_native<<<threads, blocks>>>(dev_ibead_lj_non_nat, dev_jbead_lj_non_nat, dev_itype_lj_non_nat, dev_jtype_lj_non_nat, dev_unc_pos,
                                                dev_lj_non_nat_pdb_dist, dev_value, boxl, N);

  // Free memory used for calculating binary values
  cudaFree(dev_unc_pos);
  

  // Run Scan on array
  hier_ks_scan(dev_value, dev_output, N, 0);

  int arrSize_rep;
  int endVal_rep;
  cudaMemcpy(&arrSize_rep, &dev_output[N-1], sizeof(int), cudaMemcpyDeviceToHost); 
  cudaMemcpy(&endVal_rep, &dev_value[N-1], sizeof(int), cudaMemcpyDeviceToHost); 

  // Increment arrSize by 1 if needed
  if(endVal_rep){
      arrSize_rep++;
  }

  // Copy elements
  int threads = (int)min(N, SECTION_SIZE);
  int blocks = (int)ceil(1.0*N/SECTION_SIZE);

  // Kernel to copy elements from dev_index to dev_output if their corresponding dev_value is 1
  copyElements<<<blocks, threads>>>(dev_ibead_lj_non_nat, dev_value, dev_output, dev_ibead_neighbor_list_rep, N);
  copyElements<<<blocks, threads>>>(dev_jbead_lj_non_nat, dev_value, dev_output, dev_jbead_neighbor_list_rep, N);
  copyElements<<<blocks, threads>>>(dev_jbead_lj_non_nat, dev_value, dev_output, dev_itype_neighbor_list_rep, N);
  copyElements<<<blocks, threads>>>(dev_ibead_lj_non_nat, dev_value, dev_output, dev_jtype_neighbor_list_rep, N);

  cudaFree(dev_value);
  cudaFree(dev_output);

  // Copy from device to host
  cudaMemcpy(ibead_neighbor_list_rep, dev_ibead_neighbor_list_rep, arrSize_rep*sizeof(int), cudaMemcpyDeviceToHost);
  cudaMemcpy(jbead_neighbor_list_rep, dev_jbead_neighbor_list_rep, arrSize_rep*sizeof(int), cudaMemcpyDeviceToHost);
  cudaMemcpy(itype_neighbor_list_rep, dev_itype_neighbor_list_rep, arrSize_rep*sizeof(int), cudaMemcpyDeviceToHost);
  cudaMemcpy(jtype_neighbor_list_rep, dev_jtype_neighbor_list_rep, arrSize_rep*sizeof(int), cudaMemcpyDeviceToHost);

  // Free the rest of the memory
  cudaFree(dev_ibead_lj_non_nat);
  cudaFree(dev_jbead_lj_non_nat);
  cudaFree(dev_itype_lj_non_nat);
  cudaFree(dev_jtype_lj_non_nat);
  cudaFree(dev_lj_non_nat_pdb_dist);
}

void calculate_array_native(int *&dev_ibead_lj_nat, int *&dev_jbead_lj_nat, int *&dev_itype_lj_nat, int *&dev_jtype_lj_nat, int *&dev_unc_pos, int *&dev_lj_nat_pdb_dist, 
                            int *&dev_value, int boxl, int N){
  int i = blockIdx.x * blockDim.x + threadIdx.x;
  if(i < N){
    double dx, dy, dz;
    double d2;
    int ibead, jbead, itype, jtype;
    double rcut, rcut2;

    // record sigma for ibead and jbead
    ibead = dev_ibead_lj_nat[i];
    jbead = dev_jbead_lj_nat[i];

    // record type of bead for ibead and jbead
    itype = dev_itype_lj_nat[i];
    jtype = dev_jtype_lj_nat[i];
    
    // calculate distance in x, y, and z for ibead and jbead
    dx = dev_unc_pos[jbead].x - dev_unc_pos[ibead].x;
    dy = dev_unc_pos[jbead].y - dev_unc_pos[ibead].y;
    dz = dev_unc_pos[jbead].z - dev_unc_pos[ibead].z;

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
    rcut = 3.2*sigma_rep[itype][jtype];

    // square cutoff distance, since sqrt(d2) is computationally expensive
    rcut2 = rcut*rcut;

    if(d2 < rcut2){
      dev_value[i] = 1;
    }else{
      dev_value[i] = 0;
    }
  }else if(i == 0){
    dev_value[0] = 0;
  }
}

void calculate_array_non_native(int *&dev_ibead_lj_non_nat, int *&dev_jbead_lj_non_nat, int *&dev_itype_lj_non_nat, int *&dev_jtype_lj_non_nat, 
                                int *&dev_unc_pos, int *&dev_lj_non_nat_pdb_dist, int *&dev_value, int boxl, int N){
  int i = blockIdx.x * blockDim.x + threadIdx.x+1;
  if(i < N){
    double dx, dy, dz;
    double d2;
    int ibead, jbead, itype, jtype;
    double rcut, rcut2;

    // record sigma for ibead and jbead
    ibead = dev_ibead_lj_nat[i];
    jbead = dev_jbead_lj_nat[i];

    // record type of bead for ibead and jbead
    itype = dev_itype_lj_nat[i];
    jtype = dev_jtype_lj_nat[i];
    
    // calculate distance in x, y, and z for ibead and jbead
    dx = dev_unc_pos[jbead].x - dev_unc_pos[ibead].x;
    dy = dev_unc_pos[jbead].y - dev_unc_pos[ibead].y;
    dz = dev_unc_pos[jbead].z - dev_unc_pos[ibead].z;

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
    rcut = 3.2*dev_lj_nat_pdb_dist[i];

    // square cutoff distance, since sqrt(d2) is computationally expensive
    rcut2 = rcut*rcut;

    if(d2 < rcut2){
      dev_value[i] = 1;
    }else{
      dev_value[i] = 0;
    }
  }
}

/*
 * Function: compact
 * -----------------
 *  Finds points in index with a 1 in value and stores them
 *
 *  index: array of indices to check
 *  value: binary value indicating if the corresponding index value is true (1) or false (0)
 *  N: number of elements in index and value
 *  result: pointer where compacted array is stored
 *
 *  Returns: arrSize, the size of the compacted array
 *           Note: result is modified in-place
 */

int compact(int *&dev_index, int *&dev_value, int N, int *&dev_result){
    // Perform hierarchical Kogge-Stone scan on dev_value array and store result in dev_output
    hier_ks_scan(dev_value, dev_output, N, 0);

    // Copy size of compacted array from device to host and store in arrSize
    /* 
     * TODO: If the entire array has 1 as the value, an exclusive scan will have N-1 as the last value in the array.
     * However, allocating an array with N-1 entries will not store all N values from the index array.
     * Change code to determine when we need to increment arrSize and when we don't.
     * Options include:
     *  1) Changing the hierarchical scan kernel to determine if the final value in the value array is 1
     *  2) Checking to see if the final value is 1 in the value array
     * Option 2 was selected, but please double-check this approach
     */ 
    int arrSize_att;
    int endVal_att;
    cudaMemcpy(&arrSize_att, &dev_output[N-1], sizeof(int), cudaMemcpyDeviceToHost); 
    cudaMemcpy(&endVal_att, &dev_value[N-1], sizeof(int), cudaMemcpyDeviceToHost); 

    // Increment arrSize by 1 if needed
    if(endVal_att){
        arrSize_att++;
    }

    // Declare and allocate dev_result array to store compacted indices on device (on GPU)
    int *dev_result;
    cudaMalloc((void**)&dev_result, arrSize*sizeof(int));

    // Declare and allocate dev_index to store indecies (on GPU)
    int *dev_index;
    cudaMalloc((void**)&dev_index, size);

    // Copy indices from host to device
    cudaMemcpy(dev_index, index, size, cudaMemcpyHostToDevice);

    /* Calculate number of threads and blocks to use for copying
     * If N < SECTION_SIZE (max # of threads per block), use N threads per block. Else, use SECTION_SIZE threads per block
     * Divides number of elements in array by SECTION_SIZE and rounds up, ensuring it uses the minimum number of blocks required
     */
    int threads = (int)min(N, SECTION_SIZE);
    int blocks = (int)ceil(1.0*N/SECTION_SIZE);

    // Kernel to copy elements from dev_index to dev_output if their corresponding dev_value is 1
    copyElements<<<blocks, threads>>>(dev_index, dev_value, dev_output, dev_result, N);
    
    // Sync device to ensure GPU computation is finished before proceeding
    cudaDeviceSynchronize();

    // Allocate result array on host
    result = (int *)malloc(arrSize*sizeof(int));

    // Copy dev_result (compacted array of indices in GPU) to result array on host
    cudaMemcpy(result, dev_result, arrSize*sizeof(int), cudaMemcpyDeviceToHost); 
    
    // Free device memory
    cudaFree(dev_result); 
    cudaFree(dev_index);
    cudaFree(dev_value);
    cudaFree(dev_output);

    return arrSize;
}

/*
 * Function: copyElements
 * -----------------
 *  Copys values marked true (1) from index array to result array
 *
 *  dev_index: array of indices to check (on GPU)
 *  dev_value: binary value indicating if the corresponding dev_index value is true (1) or false (0) (on GPU)
 *  N: number of elements in dev_index and dev_value
 *  dev_result: pointer where compacted array is stored (on GPU)
 */

__global__ void copyElements(int *dev_index, int *dev_value, int *dev_output, int *dev_result, int N){
    int i = blockIdx.x * blockDim.x + threadIdx.x+1;
    if(dev_value[i] && i < N){
        dev_result[dev_output[i]-1] = dev_index[i];
    }
    return;
}

/*
 * Function: hier_ks_scan
 * -----------------
 *  
 *
 *  dev_index: array of indices to check (on GPU)
 *  dev_value: binary value indicating if the corresponding dev_index value is true (1) or false (0) (on GPU)
 *  N: number of elements in dev_index and dev_value
 *  dev_result: pointer where compacted array is stored (on GPU)
 */

void hier_ks_scan(int *dev_X, int *dev_Y, int N, int re){
    if(N <= SECTION_SIZE){
        ksScanInc<<<1, N>>>(dev_X, dev_Y, N);

        cudaDeviceSynchronize();

        return;
    }else{
        int threads = (int)min(N, SECTION_SIZE);
        int blocks = (int)ceil(1.0*N/SECTION_SIZE);

        int *dev_S;
        cudaMalloc((void**)&dev_S, (int)ceil(1.0*N/SECTION_SIZE) * sizeof(int));
        
        ksScanAuxInc<<<blocks, threads>>>(dev_X, dev_Y, N, dev_S);
        cudaDeviceSynchronize();

        hier_ks_scan(dev_S, dev_S, (int)ceil(1.0*N/SECTION_SIZE), 1);
        cudaDeviceSynchronize();
        
        sumIt<<<blocks, threads>>>(dev_Y, dev_S, N);
        cudaDeviceSynchronize();

        cudaFree(dev_S);

        return;
    }
}

__global__ void ksScanAuxExc (int *X, int *Y, int InputSize, int *S) {
    int val;
    
    __shared__ int XY[SECTION_SIZE];
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if(i < InputSize && threadIdx.x != 0){
        XY[threadIdx.x] = X[i-1];
    }else{
        XY[threadIdx.x] = 0;
    }

    for(unsigned int stride = 1; stride < blockDim.x; stride *=2){
        __syncthreads();
        if(threadIdx.x >= stride){
            val = XY[threadIdx.x - stride];
        }
        __syncthreads();
        if(threadIdx.x >= stride){
            XY[threadIdx.x] += val;
        }
    }

    __syncthreads();
    if(i < InputSize){
        Y[i] = XY[threadIdx.x];
    }
    
    __syncthreads();
    if(threadIdx.x == 0){
        S[blockIdx.x] = XY[SECTION_SIZE-1];
    }
}

__global__ void ksScanAuxInc (int *X, int *Y, int InputSize, int *S) {
    int val;
    
    __shared__ int XY[SECTION_SIZE];
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if(i < InputSize){
        XY[threadIdx.x] = X[i];
    }

    for(unsigned int stride = 1; stride < blockDim.x; stride *=2){
        __syncthreads();
        if(threadIdx.x >= stride){
            val = XY[threadIdx.x - stride];
        }
        __syncthreads();
        if(threadIdx.x >= stride){
            XY[threadIdx.x] += val;
        }
    }

    __syncthreads();
    if(i < InputSize){
        Y[i] = XY[threadIdx.x];
    }
    
    __syncthreads();
    if(threadIdx.x == 0){
        S[blockIdx.x] = XY[SECTION_SIZE-1];
    }
}

__global__ void ksScanExc (int *X, int *Y, int InputSize) {
    int val;
    
    __shared__ int XY[SECTION_SIZE];
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    
    if(i < InputSize && threadIdx.x != 0){
        XY[threadIdx.x] = X[i-1];
    }else{
        XY[threadIdx.x] = 0;
    }

    for(unsigned int stride = 1; stride < blockDim.x; stride *=2){
        __syncthreads();
        if(threadIdx.x >= stride){
            val = XY[threadIdx.x - stride];
        }
        __syncthreads();
        if(threadIdx.x >= stride){
            XY[threadIdx.x] += val;
        }
    }

    __syncthreads();
    if(i < InputSize){
        Y[i] = XY[threadIdx.x];
    }
}

__global__ void ksScanInc (int *X, int *Y, int InputSize) {
    int val;
    
    __shared__ int XY[SECTION_SIZE];
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    
    if(i < InputSize){
        XY[threadIdx.x] = X[i];
    }

    for(unsigned int stride = 1; stride < blockDim.x; stride *=2){
        __syncthreads();
        if(threadIdx.x >= stride){
            val = XY[threadIdx.x - stride];
        }
        __syncthreads();
        if(threadIdx.x >= stride){
            XY[threadIdx.x] += val;
        }
    }

    __syncthreads();
    if(i < InputSize){
        Y[i] = XY[threadIdx.x];
    }
}

__global__ void sumIt (int *Y, int *S, int InputSize) {
    if(blockIdx.x > 0){
        int i = blockIdx.x * blockDim.x + threadIdx.x;
        if(i < InputSize){
            Y[i] += S[blockIdx.x-1];
        }
    }
}

__global__ void init_kernel(){
 return;
}