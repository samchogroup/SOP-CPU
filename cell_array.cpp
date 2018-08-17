#include <math.h>
#include <cstdlib>
#include <cstring>
#include <cstdio>
#include "cell_array.h"
#include "global.h"

void update_cell_array() {
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

	//printf("X range: %f, Y range: %f, Z range: %f\n", maxX - minX, maxY - minY, maxZ - minZ);
	
	//printf("X range: %f, Y range: %f, Z range: %f\n", maxX - minX, maxY - minY, maxZ - minZ);
	//printf("X cell: %f, Y cell: %f, Z cell: %f\n", lcellX, lcellY, lcellZ);

	//int numCells = (ncell)*(ncell)*(ncell);
	//int * beadLinks = new int[nbead];
	//int * cells = new int[numCells];
	
	memset(cells, -1, numCells * sizeof(int));
	memset(beadLinks, -1, nbead * sizeof(int));
	
	//int offset = floor((ncell-1)/2.0);
	//printf("Box length, num cells, length of cell = %f, %f, %f\n", boxl,ncell,lcell);
	for (int i = 1; i <= nbead; i++ ) {
		
		imcx = floor((pos[i].x - minX) / lcell) ;
		imcy = floor((pos[i].y - minY) / lcell) ;
		imcz = floor((pos[i].z - minZ) / lcell) ;
		
		index = (imcz) * (ncell)*(ncell) + imcy * (ncell) + imcx;
		//printf("Bead %d has cell dimensions (%f, %f, %f).\n", i, pos[i].x, pos[i].y, pos[i].z);
		//printf("Bead %d has cell indices (%d, %d, %d) and index %d.\n", i, imcx, imcy, imcz, index);

		beadLinks[i] = cells[index];
		cells[index] = i;

	}
	//exit(0);
	nil_att = 0;
 	nil_rep = 0;
 	/*
 	for (int i = 1; i<= nbead; i++) {
 		printf("Bead %d has cell dimensions (%f, %f, %f).\n", i, pos[i].x, pos[i].y, pos[i].z);
 	}
	//exit(0);
	*/
	int tempx, tempy, tempz, cellN;
	cellN = int(ncell);
	for (int i = 1; i<= nbead - 3; i++) {

		imcx = floor((pos[i].x - minX) / lcell) ;
		imcy = floor((pos[i].y - minY) / lcell) ;
		imcz = floor((pos[i].z - minZ) / lcell) ;
    	//index = (imcz) * (ncell)*(ncell) + imcy * (ncell) + imcx;
		displacementBead = (nbead-3)*(i-1) - (i-2)*(i-1)/2;
		//printf("Bead %d has cell dimensions (%f, %f, %f).\n", i, pos[i].x, pos[i].y, pos[i].z);
		//printf("Bead %d has cell indices (%d, %d, %d) and index %d.\n", i, imcx, imcy, imcz, index);

		for(int deltaZ = -1; deltaZ <= 1; deltaZ++) {
			for(int deltaY = -1; deltaY <= 1; deltaY++) {
				for(int deltaX = -1; deltaX <= 1; deltaX++) {
					tempz = (imcz + deltaZ + cellN)%(cellN);
					tempy = (imcy + deltaY + cellN)%(cellN);
					tempx = (imcx + deltaX + cellN)%(cellN);
					
					index = tempz * ncell*ncell + tempy * ncell + tempx;
					//printf("Bead %d has cell indices (%d, %d, %d) and index %d.\n", i, tempx, tempy, tempz, index);
					j = cells[index];					

					while(j != -1 && j > i+2) {
						displacement = displacementBead +  j - 2 - i;

						itype = itype_lj_tot[displacement];
						jtype = jtype_lj_tot[displacement];

						//printf("Bead %d has neighbor %d with indices (%d, %d, %d)\n", i, j, tempx, tempy, tempz);
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
	//printf("Leaving cell_array! %d; %d, %d\n", nil_att, nil_rep, ibead_pair_list_rep[1]);
	//exit(0);
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
