#include <math.h>
#include <cstdlib>
#include <cstring>
#include <cstdio>
#include "cell_array.h"
#include "global.h"
/*
Function used to sort beads into cells and then place interactions between cells into the pair list using this cell array structure.
Interactions are taken from the total lj list and placed into the pair list for force evaluation.
*/
void update_cell_array() {
	int imcx, imcy, imcz, index;
	int itype, jtype, j, displacementBead, displacement;
	
	double minX, minY, minZ;

	// Get x,y,z, position of first bead in position array.
	minX = pos[1].x;
	minY = pos[1].y;
	minZ = pos[1].z;

	// Loop through all beads and find minimum positions.
	// Minimum position is needed when constructing cell array so that we work with positive coordinates.
	for(int i = 2; i <= nbead; i++) {
		if (pos[i].x < minX) {
			minX = pos[i].x;
		}

		if (pos[i].y < minY) {
			minY = pos[i].y;
		}

		if (pos[i].z < minZ) {
			minZ = pos[i].z;
		}
	}

	
	// Set the cells and beadLinks array to have -1 entries in each location. -1 means that no bead has yet been added to that cell
	memset(cells, -1, numCells * sizeof(int));
	memset(beadLinks, -1, nbead * sizeof(int));
	
	// Loop through all beads (no 0th bead) and place them in a cell
	for (int i = 1; i <= nbead; i++ ) {
		
		// Calculate the index of the cell the bead belongs to. One index for each direction.
		imcx = floor((pos[i].x - minX) / lcell) ;
		imcy = floor((pos[i].y - minY) / lcell) ;
		imcz = floor((pos[i].z - minZ) / lcell) ;
		
		// Combine three indices into one index that is specific to the cell.
		index = (imcz) * (ncell)*(ncell) + imcy * (ncell) + imcx;

		// Have current bead point to the last bead in the cell. Then make current bead the last one in the cell specificed by index.
		beadLinks[i] = cells[index];
		cells[index] = i;

	}
	// We now have constructed the cell array structure. Each bead has been placed in a cell and points to the next bead in that cell.
	// Essentially we have a singly linked list implemented using an array.
	
	// Set the number of attractive and repulsive interactions placed in pair list to 0.
	nil_att = 0;
 	nil_rep = 0;
	
	int tempx, tempy, tempz, cellN;
	cellN = int(ncell);

	// Bead interactions are only calculated for beads i and those greater than i + 2.
	// i and i+1 interactions are taken care of by bonds.
	// i and i+2 are taken care of by dihedral angles.
	for (int i = 1; i<= nbead - 3; i++) {
		// Calcualte the cell indices of the bead current being consider.
		imcx = floor((pos[i].x - minX) / lcell) ;
		imcy = floor((pos[i].y - minY) / lcell) ;
		imcz = floor((pos[i].z - minZ) / lcell) ;
		
		// Calculate the index of the first interaction beginning with this bead in the list of all lj interactions.
		displacementBead = (nbead-3)*(i-1) - (i-2)*(i-1)/2;

		// Loop through the surrounding 27 cells and pull out the interactions between beads in these cells and the current bead.
		for(int deltaZ = -1; deltaZ <= 1; deltaZ++) {
			for(int deltaY = -1; deltaY <= 1; deltaY++) {
				for(int deltaX = -1; deltaX <= 1; deltaX++) {
					// The modding by cellN handles periodic boundary conditions.
					tempz = (imcz + deltaZ + cellN)%(cellN);
					tempy = (imcy + deltaY + cellN)%(cellN);
					tempx = (imcx + deltaX + cellN)%(cellN);
					
					// Calcualte the single index for the neighboring cell.
					index = tempz * ncell*ncell + tempy * ncell + tempx;
					// Get first bead in this neighboring cell
					j = cells[index];					

					// Loop through all beads in this cell, only looking at those with a higher bead number.
					// -1 indicates the end of a cell
					// Bead numbers will decrease as we move in a cell, so we stop when we incounter bead i+2.
					while(j != -1 && j > i+2) {
						displacement = displacementBead +  j - 2 - i;

						itype = itype_lj_tot[displacement];
						jtype = jtype_lj_tot[displacement];

						if (lg_tot_pdb_dist[displacement] < rcut_nat[itype][jtype]) {
							// Place interaction in attractive Pair List
							nil_att++;
							ibead_pair_list_att[nil_att] = ibead_lj_tot[displacement];
							jbead_pair_list_att[nil_att] = jbead_lj_tot[displacement];
							itype_pair_list_att[nil_att] = itype;
							jtype_pair_list_att[nil_att] = jtype;

							pl_lj_nat_pdb_dist[nil_att] = lg_tot_pdb_dist[displacement];
							pl_lj_nat_pdb_dist2[nil_att] = lg_tot_pdb_dist2[displacement];
							pl_lj_nat_pdb_dist6[nil_att] = lg_tot_pdb_dist6[displacement];
							pl_lj_nat_pdb_dist12[nil_att] = lg_tot_pdb_dist12[displacement];
						}
						else {
							// Place interaction in repulsive Pair List
							nil_rep++;
							
							ibead_pair_list_rep[nil_rep] = ibead_lj_tot[displacement];
							jbead_pair_list_rep[nil_rep] = jbead_lj_tot[displacement];
							itype_pair_list_rep[nil_rep] = itype;
							jtype_pair_list_rep[nil_rep] = jtype;
						}
						// Get next bead in that cell.
						j = beadLinks[j];
					}					
				}
			}		
		}
		// Done looping through all 27 surrounding cells.
	}
	// Done looping through all beads.
	// All interactions to be considered for energy and position updates have been placed in Pair List.
}


/* Same as above function, but rather than put beads into pair list, put into neighbor list. */
void update_cell_array_hybrid() {
	int imcx, imcy, imcz, index;
	int itype, jtype, j, displacementBead, displacement;
	
	double minX, minY, minZ;
	
	minX = pos[1].x;
	minY = pos[1].y;
	minZ = pos[1].z;

	for(int i = 2; i <= nbead; i++) {
		if (pos[i].x < minX) {
			minX = pos[i].x;
		}

		if (pos[i].y < minY) {
			minY = pos[i].y;
		}

		if (pos[i].z < minZ) {
			minZ = pos[i].z;
		}
	}
	
	memset(cells, -1, numCells * sizeof(int));
	memset(beadLinks, -1, nbead * sizeof(int));
	
	for (int i = 1; i <= nbead; i++ ) {
		
		imcx = floor((pos[i].x - minX) / lcell) ;
		imcy = floor((pos[i].y - minY) / lcell) ;
		imcz = floor((pos[i].z - minZ) / lcell) ;
		
		index = (imcz) * (ncell)*(ncell) + imcy * (ncell) + imcx;
		
		beadLinks[i] = cells[index];
		cells[index] = i;

	}

	nnl_att = 0;
 	nnl_rep = 0;
 	

	int tempx, tempy, tempz, cellN;
	cellN = int(ncell);
	for (int i = 1; i<= nbead - 3; i++) {

		imcx = floor((pos[i].x - minX) / lcell) ;
		imcy = floor((pos[i].y - minY) / lcell) ;
		imcz = floor((pos[i].z - minZ) / lcell) ;

		displacementBead = (nbead-3)*(i-1) - (i-2)*(i-1)/2;
		
		for(int deltaZ = -1; deltaZ <= 1; deltaZ++) {
			for(int deltaY = -1; deltaY <= 1; deltaY++) {
				for(int deltaX = -1; deltaX <= 1; deltaX++) {
					tempz = (imcz + deltaZ + cellN)%(cellN);
					tempy = (imcy + deltaY + cellN)%(cellN);
					tempx = (imcx + deltaX + cellN)%(cellN);
					
					index = tempz * ncell*ncell + tempy * ncell + tempx;

					j = cells[index];					

					while(j != -1 && j > i+2) {
						displacement = displacementBead +  j - 2 - i;

						itype = itype_lj_tot[displacement];
						jtype = jtype_lj_tot[displacement];

						if (lg_tot_pdb_dist[displacement] < rcut_nat[itype][jtype]) {
							nnl_att++;

							ibead_neighbor_list_att[nnl_att] = ibead_lj_tot[displacement];
							jbead_neighbor_list_att[nnl_att] = jbead_lj_tot[displacement];
							itype_neighbor_list_att[nnl_att] = itype;
							jtype_neighbor_list_att[nnl_att] = jtype;

							nl_lj_nat_pdb_dist[nnl_att] = lg_tot_pdb_dist[displacement];
							nl_lj_nat_pdb_dist2[nnl_att] = lg_tot_pdb_dist2[displacement];
							nl_lj_nat_pdb_dist6[nnl_att] = lg_tot_pdb_dist6[displacement];
							nl_lj_nat_pdb_dist12[nnl_att] = lg_tot_pdb_dist12[displacement];

						}
						else {
							nnl_rep++;
							
							ibead_neighbor_list_rep[nnl_rep] = ibead_lj_tot[displacement];
							jbead_neighbor_list_rep[nnl_rep] = jbead_lj_tot[displacement];
							itype_neighbor_list_rep[nnl_rep] = itype;
							jtype_neighbor_list_rep[nnl_rep] = jtype;
							
						}
						j = beadLinks[j];
					}					
				}
			}		
		}
	}
}

void update_cell_array2() {
	int imcx, imcy, imcz, index;
	int itype, jtype, j, displacementBead, displacement;
	
	double maxX,maxY,maxZ, minX, minY, minZ, lcellX, lcellY, lcellZ;
	
	maxX = minX = unc_pos[1].x;
	maxY = minY = unc_pos[1].y;
	maxZ = minZ = unc_pos[1].z;

	for(int i = 2; i <= nbead; i++) {
		if (unc_pos[i].x > maxX) {
			maxX = unc_pos[i].x;
		}
		else if (unc_pos[i].x < minX) {
			minX = unc_pos[i].x;
		}

		if (unc_pos[i].y > maxY) {
			maxY = unc_pos[i].y;
		}
		else if (unc_pos[i].y < minY) {
			minY = unc_pos[i].y;
		}

		if (unc_pos[i].z > maxZ) {
			maxZ = unc_pos[i].z;
		}
		else if (unc_pos[i].z < minZ) {
			minZ = unc_pos[i].z;
		}
	}

	lcellX = (maxX - minX) / double(ncell);
	lcellY = (maxY - minY) / double(ncell);
	lcellZ = (maxZ - minZ) / double(ncell);
	
	//printf("X range: %f, Y range: %f, Z range: %f\n", maxX - minX, maxY - minY, maxZ - minZ);
	//printf("X cell: %f, Y cell: %f, Z cell: %f\n", lcellX, lcellY, lcellZ);

	int numCellsNEW = (ncell+2)*(ncell+2)*(ncell+2);
	//int numCellsNEW = (ncell)*(ncell)*(ncell);
	int * beadLinks = new int[nbead];
	int * cells = new int[numCellsNEW];
	memset(cells, -1, numCellsNEW * sizeof(int));
	memset(beadLinks, -1, nbead * sizeof(int));
	
	int offset = floor((ncell-1)/2.0);
	//printf("Box length, num cells, length of cell = %f, %f, %f\n", boxl,ncell,lcell);
	for (int i = 1; i <= nbead; i++ ) {
		
		imcx = (int(unc_pos[i].x - minX - 0.1) / lcellX)+1;
	   	imcy = (int(unc_pos[i].y - minY - 0.1) / lcellY)+1;
	    imcz = (int(unc_pos[i].z - minZ - 0.1) / lcellZ)+1;
	    
	    /*
	    imcx = ((unc_pos[i].x - minX) / lcellX);
	   	imcy = ((unc_pos[i].y - minY) / lcellY);
	    imcz = ((unc_pos[i].z - minZ) / lcellZ);
		*/
		index = (imcz) * (ncell+2)*(ncell+2) + imcy * (ncell+2) + imcx;
		//index = (imcz) * (ncell)*(ncell) + imcy * (ncell) + imcx;
		//printf("Bead %d has cell dimensions (%f, %f, %f).\n", i, unc_pos[i].x, unc_pos[i].y, unc_pos[i].z);
		// printf("Bead %d has cell indices (%d, %d, %d) and index %d.\n", i, imcx, imcy, imcz, index);
		
			beadLinks[i] = cells[index];
			cells[index] = i;
	
	}
	// exit(0);
	nil_att = 0;
 	nil_rep = 0;
 	/*
 	for (int i = 1; i<= nbead; i++) {
 		printf("Bead %d has cell dimensions (%f, %f, %f).\n", i, unc_pos[i].x, unc_pos[i].y, unc_pos[i].z);
 	}
	exit(0);
	*/
	int tempx, tempy, tempz,cellN;
	cellN = int(ncell);
	for (int i = 1; i<= nbead - 3; i++) {
		
		
		imcx = (int(unc_pos[i].x - minX - 0.1) / lcellX) + 1;
	   	imcy = (int(unc_pos[i].y - minY - 0.1) / lcellY) + 1;
	   	imcz = (int(unc_pos[i].z - minZ - 0.1) / lcellZ) + 1;
	   	
	   	
	   	/*
	   	imcx = (int(unc_pos[i].x - minX - 0.1) / lcellX) ;
	   	imcy = (int(unc_pos[i].y - minY - 0.1) / lcellY) ;
	   	imcz = (int(unc_pos[i].z - minZ - 0.1) / lcellZ) ;
		*/
		// ibead = ibead_lj_nat[i];
    		// itype = itype_lj_nat[i];
    	
		displacementBead = (nbead-3)*(i-1) - (i-2)*(i-1)/2;
		//printf("Bead %d has cell dimensions (%f, %f, %f).\n", i, unc_pos[i].x, unc_pos[i].y, unc_pos[i].z);
		//printf("Bead %d has cell indices (%d, %d, %d)\n", i, imcx, imcy, imcz);

		for(int deltaZ = -1; deltaZ <= 1; deltaZ++) {
			for(int deltaY = -1; deltaY <= 1; deltaY++) {
				for(int deltaX = -1; deltaX <= 1; deltaX++) {
					//tempz = (imcz + deltaZ + cellN)%(cellN);
					//tempy = (imcy + deltaY + cellN)%(cellN);
					//tempx = (imcx + deltaX + cellN)%(cellN);
					
					// index = ((imcz + deltaZ)%int(ncell)) * ncell*ncell + ((imcy + deltaY)%int(ncell)) * ncell + (imcx + deltaX)%int(ncell);
					//index = tempz * ncell*ncell + tempy * ncell + tempx;
					// printf("Bead %d has cell indices (%d, %d, %d) and index %d.\n", i, tempx, tempy, tempz, index);
					index = (imcz + deltaZ) * (ncell+2)*(ncell+2) + (imcy + deltaY) * (ncell+2) + (imcx + deltaX);
					j = cells[index];					

					while(j != -1 && j > i+2) {
						//jbead = jbead_lj_nat[i];
						//jtype = jtype_lj_nat[i];
						displacement = displacementBead +  j - 2 - i;

						// ibead = ibead_lj_nat[displacement];
    						// jbead = jbead_lj_nat[displacement];
    						itype = itype_lj_tot[displacement];
    						jtype = jtype_lj_tot[displacement];

						//printf("Bead %d has neighbor %d with indices (%d, %d, %d)\n", i, j, imcx + deltaX, imcy + deltaY, imcz + deltaZ);
						//printf("Displacement: %d; VDW: %f\n", displacement, lg_tot_pdb_dist[displacement]);
						//printf("i type: %d; j type: %d\n", itype, jtype);
						//printf("First entry in repulsive ibead array: %d %d\n", nil_rep, ibead_pair_list_rep[1]);
						if (lg_tot_pdb_dist[displacement] < rcut_nat[itype][jtype]) {
							//printf("Here Attractive\n");
							nil_att++;
							ibead_pair_list_att[nil_att] = ibead_lj_tot[displacement];
							jbead_pair_list_att[nil_att] = jbead_lj_tot[displacement];
							itype_pair_list_att[nil_att] = itype;
							jtype_pair_list_att[nil_att] = jtype;

							pl_lj_nat_pdb_dist[nil_att] = lg_tot_pdb_dist[displacement];
							pl_lj_nat_pdb_dist2[nil_att] = lg_tot_pdb_dist2[displacement];
							pl_lj_nat_pdb_dist6[nil_att] = lg_tot_pdb_dist6[displacement];
							pl_lj_nat_pdb_dist12[nil_att] = lg_tot_pdb_dist12[displacement];
							//printf("Here Attractive2\n");
						}
						else {
							//printf("Here Repulsive\n");
							nil_rep++;
							
							ibead_pair_list_rep[nil_rep] = ibead_lj_tot[displacement];
							//printf("Displacement: %d; %d\n", ibead_lj_tot[displacement], ibead_pair_list_rep[nil_rep]);
							//exit(0);
							jbead_pair_list_rep[nil_rep] = jbead_lj_tot[displacement];
							itype_pair_list_rep[nil_rep] = itype;
							jtype_pair_list_rep[nil_rep] = jtype;
							//printf("%d, %d, %d, %d\n", ibead_pair_list_rep[nil_rep], jbead_pair_list_rep[nil_rep], itype_pair_list_rep[nil_rep], jtype_pair_list_rep[nil_rep]);
							//exit(0);
							//printf("Here Repulsive2\n");
						}
						j = beadLinks[j];
						//printf("New bead: %d\n", j);
						//printf("First entry in repulsive ibead array: %d %d\n", nil_rep, ibead_pair_list_rep[1]);
					}					
				}
				//printf("First entry in repulsive ibead array: %d %d\n", nil_rep, ibead_pair_list_rep[1]);	
			}		
		}
	}
	//ibead_pair_list_rep[1] = 1;
	// printf("Leaving cell_array! %d; %d, %d\n", nil_att, nil_rep, ibead_pair_list_rep[1]);
	//exit(0);
}
