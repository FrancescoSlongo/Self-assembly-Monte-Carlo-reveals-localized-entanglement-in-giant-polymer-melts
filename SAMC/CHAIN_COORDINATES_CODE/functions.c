///////////////////////////////////////////////////////////
// efornasa, February 2026                               //
// Implementation of the functions defined in the header //
// file "functions.h" for extracting chain coordinates   //
// from configurations sampled with SAMC.                //
///////////////////////////////////////////////////////////
// construct_lattice.................................... //
// construct_lattice_from_input_file.................... //
// free_lattice......................................... //
// extract_linear_chain................................. //
// extract_ring......................................... //
// find_next_site....................................... //
// append_element_to_int_array.......................... //
// free_chain........................................... //
// construct_configuration_chains....................... //
// free_configuration................................... //
///////////////////////////////////////////////////////////

// functions.c

#include <stdio.h>
#include <stdbool.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "functions.h"


//////////////////////////////////////////////////////////////////////////////
// HYPER_CUBIC_LATTICE CONSTRUCTOR
HC_Lattice *construct_lattice(int dimension, int size) {
  // Allocate memory for the HC_Lattice structure
  HC_Lattice *lattice = (HC_Lattice *)malloc(sizeof(HC_Lattice));
  if (!lattice) {
    fprintf(stderr, "Memory allocation failed for HC_Lattice.\n");
    return NULL;
  }

  // Initialize the dimension and size
  lattice->dimension = dimension;
  lattice->size = size;

  // Set the number of sites and edges
  lattice->n_sites = pow(size, dimension);
  lattice->n_edges = dimension * lattice->n_sites;

  // ALLOCATE ARRAYS' MEMORY
	int i, j;

  // Allocate memory for the site array
  lattice->site = (int **)malloc(lattice->n_sites * sizeof(int *));
  if (!lattice->site) {
    fprintf(stderr, "Memory allocation failed for site array.\n");
    free(lattice);
    return NULL;
  }
  for (i = 0; i < lattice->n_sites; i++) {
    lattice->site[i] = (int *)malloc(dimension * sizeof(int));
    if (!lattice->site[i]) {
      fprintf(stderr, "Memory allocation failed for site array row %d.\n", i);
      for (j = 0; j < i; j++) free(lattice->site[j]);
      free(lattice->site);
      free(lattice);
      return NULL;
    }
  }

	// Allocate memory for nearest-neighbor array
	lattice->nearest_neighbor = (int **)malloc(lattice->n_sites * sizeof(int *));
  if (!lattice->nearest_neighbor) {
    fprintf(stderr, "Memory allocation failed for site array.\n");
		for (i = 0; i < lattice->n_sites; i++) free(lattice->site[i]);
		free(lattice->site);
    free(lattice);
    return NULL;
	}
	for (i = 0; i < lattice->n_sites; i++) {
		lattice->nearest_neighbor[i] = (int *)malloc(2 * dimension * sizeof(int));
		if (lattice->nearest_neighbor[i] == NULL) {
	    fprintf(stderr, "Memory allocation failed for nearest-neighbor array row %d.\n", i);
			for (j = 0; j < lattice->n_sites; j++) free(lattice->site[j]);
			free(lattice->site);
			for (j = 0; j < i; j++) free(lattice->nearest_neighbor[j]);
			free(lattice->nearest_neighbor);
    	free(lattice);
    	return NULL;
	  }
	}

  // Allocate memory for edge array
  lattice->edge = (int **)malloc(lattice->n_edges * sizeof(int *));
  if (!lattice->edge) {
    fprintf(stderr, "Memory allocation failed for edge array.\n");
    for (i = 0; i < lattice->n_sites; i++) {
			free(lattice->site[i]);
			free(lattice->nearest_neighbor[i]);
		}
    free(lattice->site);
		free(lattice->nearest_neighbor);
    free(lattice);
    return NULL;
  }
  for (i = 0; i < lattice->n_edges; i++) {
    lattice->edge[i] = (int *)malloc(2 * sizeof(int));
    if (!lattice->edge[i]) {
      fprintf(stderr, "Memory allocation failed for edge array row %d.\n", i);
      for (j = 0; j < lattice->n_sites; j++) {
				free(lattice->site[j]);
				free(lattice->nearest_neighbor[j]);
			}
      free(lattice->site);
      free(lattice->nearest_neighbor);
      for (j = 0; j < i; j++) free(lattice->edge[j]);
      free(lattice->edge);
      free(lattice);
      return NULL;
    }
  }
  
  // Allocate memory for adjacent-edge array
  lattice->adjacent_edge = (int **)malloc(lattice->n_edges * sizeof(int *));
  if (!lattice->adjacent_edge) {
    fprintf(stderr, "Memory allocation failed for adjacent-edge array.\n");
    for (i = 0; i < lattice->n_sites; i++) {
			free(lattice->site[i]);
			free(lattice->nearest_neighbor[i]);
		}
    free(lattice->site);
		free(lattice->nearest_neighbor);
		for (i = 0; i < lattice->n_edges; i++) free(lattice->edge[i]);
		free(lattice->edge);
    free(lattice);
    return NULL;
  }
  for (i = 0; i < lattice->n_edges; i++) {
    lattice->adjacent_edge[i] = (int *)malloc(2 * (2 * dimension - 1) * sizeof(int));
    if (!lattice->adjacent_edge[i]) {
      fprintf(stderr, "Memory allocation failed for adjacent-dge array row %d.\n", i);
      for (j = 0; j < lattice->n_sites; j++) {
				free(lattice->site[j]);
				free(lattice->nearest_neighbor[j]);
			}
      free(lattice->site);
      free(lattice->nearest_neighbor);
			for (i = 0; i < lattice->n_edges; i++) free(lattice->edge[i]);
			free(lattice->edge);
      for (j = 0; j < i; j++) free(lattice->adjacent_edge[j]);
      free(lattice->adjacent_edge);
      free(lattice);
      return NULL;
    }
  }

  // Allocate memory for incident-edge array
  lattice->incident_edge = (int **)malloc(lattice->n_sites * sizeof(int *));
  if (!lattice->incident_edge) {
    fprintf(stderr, "Memory allocation failed for incident-edge array.\n");
    for (i = 0; i < lattice->n_sites; i++) {
      free(lattice->site[i]);
      free(lattice->nearest_neighbor[i]);
    }
    free(lattice->site);
    free(lattice->nearest_neighbor);
    for (i = 0; i < lattice->n_edges; i++) {
      free(lattice->edge[i]);
      free(lattice->adjacent_edge[i]);
    }
    free(lattice->edge);
    free(lattice->adjacent_edge);
    free(lattice);
    return NULL;
  }
  for (i = 0; i < lattice->n_sites; i++) {
    lattice->incident_edge[i] = (int *)malloc(2 * dimension * sizeof(int));
    if (!lattice->incident_edge[i]) {
      fprintf(stderr, "Memory allocation failed for incident-edge array row %d.\n", i);
      for (j = 0; j < lattice->n_sites; j++) {
				free(lattice->site[j]);
				free(lattice->nearest_neighbor[j]);
      }
      free(lattice->site);
      free(lattice->nearest_neighbor);
      for (i = 0; i < lattice->n_edges; i++) {
				free(lattice->edge[i]);
				free(lattice->adjacent_edge[i]);
      }
      free(lattice->edge);
      free(lattice->adjacent_edge);
      for (j = 0; j < i; j++) free(lattice->incident_edge[j]);
      free(lattice->incident_edge);
      free(lattice);
      return NULL;
    }
  }
  
  // POPULATE ARRAYS
  // Populate sites
  int aux;
  for (i = 0; i < lattice->n_sites; i++) {
    aux = i;
    for (j = 0; j < dimension; j++) {
      lattice->site[i][j] = aux % size;
      aux /= size;
    }
  }

	// Populate nearest-neighbors
  for (i = 0; i < lattice->n_sites; i++) {
		for (j = 0; j < dimension; j++) {
			lattice->nearest_neighbor[i][2 * j] = i - pow(size, j);
			lattice->nearest_neighbor[i][2 * j + 1] = i + pow(size, j);
			if (lattice->site[i][j] == 0) lattice->nearest_neighbor[i][2 * j] += pow(size, j + 1);
			else if (lattice->site[i][j] == (lattice->size - 1)) lattice->nearest_neighbor[i][2 * j + 1] -= pow(size, j + 1);
		}
	}	

  // Populate edges
	aux = 0;
	for (i = 0; i < lattice->n_sites; i++) {
		for (j = 0; j < 2 * dimension; j++) {
			if (lattice->nearest_neighbor[i][j] > i) {
				lattice->edge[aux][0] = i;
				lattice->edge[aux][1] = lattice->nearest_neighbor[i][j];
				aux++;
			}
		}
	}

	// Populate adjacent-edge array
	for (i = 0; i < lattice->n_edges; i++) {
		aux = 0;
		for (j = 0; j < lattice->n_edges; j++) {
			if (i != j) {
				if((lattice->edge[i][0] == lattice->edge[j][0]) || (lattice->edge[i][0] == lattice->edge[j][1]) || (lattice->edge[i][1] == lattice->edge[j][0]) || (lattice->edge[i][1] == lattice->edge[j][1])) {
					lattice->adjacent_edge[i][aux] = j;
					aux++;
				}
			}
		}
	}

  // Populate incident-edge array
  int *counter = (int *)malloc(lattice->n_sites * sizeof(int));
  for (i = 0; i < lattice->n_sites; i++) counter[i] = 0;
  for (i = 0; i < lattice->n_edges; i++) {
    aux = lattice->edge[i][0];
    lattice->incident_edge[aux][counter[aux]] = i;
    counter[aux]++;
    aux = lattice->edge[i][1];
    lattice->incident_edge[aux][counter[aux]] = i;
    counter[aux]++;
  }
  free(counter);

  return lattice;
}


//////////////////////////////////////////////////////////////////////////////
// HYPER_CUBIC_LATTICE CONSTRUCTOR FROM INPUT FILE
HC_Lattice *construct_lattice_from_input_file(const char *lattice_filename) {
  FILE *file = fopen(lattice_filename, "r");
  if (!file) {
    printf("Error: Could not open file %s\n", lattice_filename);
    return NULL;
  }
  
  // Allocate memory for the HC_Lattice structure
  HC_Lattice *lattice = (HC_Lattice *)malloc(sizeof(HC_Lattice));
  if (!lattice) {
    fprintf(stderr, "Memory allocation failed for HC_Lattice.\n");
    return NULL;
  }
  
  // Read dimension and size
  if (fscanf(file, "#dimension %d\n", &lattice->dimension) == EOF) {
    printf("Fatal error while reading input file.\n");
    exit(1);
  }
  if (fscanf(file, "#size %d\n", &lattice->size) == EOF) {
    printf("Fatal error while reading input file.\n");
    exit(1);
  }

  lattice->n_sites = pow(lattice->size, lattice->dimension);
  lattice->n_edges = lattice->dimension * lattice->n_sites;

  int i, j, tmp;

  // MEMORY ALLOCATION
  // Allocate memory for the site array
  lattice->site = (int **)malloc(lattice->n_sites * sizeof(int *));
  if (!lattice->site) {
    fprintf(stderr, "Memory allocation failed for site array.\n");
    free(lattice);
    return NULL;
  }
  for (i = 0; i < lattice->n_sites; i++) {
    lattice->site[i] = (int *)malloc(lattice->dimension * sizeof(int));
    if (!lattice->site[i]) {
      fprintf(stderr, "Memory allocation failed for site array row %d.\n", i);
      for (j = 0; j < i; j++) free(lattice->site[j]);
      free(lattice->site);
      free(lattice);
      return NULL;
    }
  }
  
  // Allocate memory for nearest-neighbor array
  lattice->nearest_neighbor = (int **)malloc(lattice->n_sites * sizeof(int *));
  if (!lattice->nearest_neighbor) {
    fprintf(stderr, "Memory allocation failed for site array.\n");
    for (i = 0; i < lattice->n_sites; i++) free(lattice->site[i]);
    free(lattice->site);
    free(lattice);
    return NULL;
  }
  for (i = 0; i < lattice->n_sites; i++) {
    lattice->nearest_neighbor[i] = (int *)malloc(2 * lattice->dimension * sizeof(int));
    if (!lattice->nearest_neighbor[i]) {
      fprintf(stderr, "Memory allocation failed for nearest-neighbor array row %d.\n", i);
      for (j = 0; j < lattice->n_sites; j++) free(lattice->site[j]);
      free(lattice->site);
      for (j = 0; j < i; j++) free(lattice->nearest_neighbor[j]);
      free(lattice->nearest_neighbor);
      free(lattice);
      return NULL;
    }
  }
  
  // Allocate memory for edge array
  lattice->edge = (int **)malloc(lattice->n_edges * sizeof(int *));
  if (!lattice->edge) {
    fprintf(stderr, "Memory allocation failed for edge array.\n");
    for (i = 0; i < lattice->n_sites; i++) {
      free(lattice->site[i]);
      free(lattice->nearest_neighbor[i]);
    }
    free(lattice->site);
    free(lattice->nearest_neighbor);
    free(lattice);
    return NULL;
  }
  for (i = 0; i < lattice->n_edges; i++) {
    lattice->edge[i] = (int *)malloc(2 * sizeof(int));
    if (!lattice->edge[i]) {
      fprintf(stderr, "Memory allocation failed for edge array row %d.\n", i);
      for (j = 0; j < lattice->n_sites; j++) {
      	free(lattice->site[j]);
      	free(lattice->nearest_neighbor[j]);
      }
      free(lattice->site);
      free(lattice->nearest_neighbor);
      for (j = 0; j < i; j++) free(lattice->edge[j]);
      free(lattice->edge);
      free(lattice);
      return NULL;
    }
  }
  
  // Allocate memory for adjacent-edge array
  lattice->adjacent_edge = (int **)malloc(lattice->n_edges * sizeof(int *));
  if (!lattice->adjacent_edge) {
    fprintf(stderr, "Memory allocation failed for adjacent-edge array.\n");
    for (i = 0; i < lattice->n_sites; i++) {
      free(lattice->site[i]);
      free(lattice->nearest_neighbor[i]);
    }
    free(lattice->site);
    free(lattice->nearest_neighbor);
    for (i = 0; i < lattice->n_edges; i++) free(lattice->edge[i]);
    free(lattice->edge);
    free(lattice);
    return NULL;
  }
  for (i = 0; i < lattice->n_edges; i++) {
    lattice->adjacent_edge[i] = (int *)malloc(2 * (2 * lattice->dimension - 1) * sizeof(int));
    if (!lattice->adjacent_edge[i]) {
      fprintf(stderr, "Memory allocation failed for adjacent-dge array row %d.\n", i);
      for (j = 0; j < lattice->n_sites; j++) {
	      free(lattice->site[j]);
	      free(lattice->nearest_neighbor[j]);
      }
      free(lattice->site);
      free(lattice->nearest_neighbor);
      for (i = 0; i < lattice->n_edges; i++) free(lattice->edge[i]);
      free(lattice->edge);
      for (j = 0; j < i; j++) free(lattice->adjacent_edge[j]);
      free(lattice->adjacent_edge);
      free(lattice);
      return NULL;
    }
  }
  
  // Allocate memory for incident-edge array
  lattice->incident_edge = (int **)malloc(lattice->n_sites * sizeof(int *));
  if (!lattice->incident_edge) {
    fprintf(stderr, "Memory allocation failed for incident-edge array.\n");
    for (i = 0; i < lattice->n_sites; i++) {
      free(lattice->site[i]);
      free(lattice->nearest_neighbor[i]);
    }
    free(lattice->site);
    free(lattice->nearest_neighbor);
    for (i = 0; i < lattice->n_edges; i++) {
      free(lattice->edge[i]);
      free(lattice->adjacent_edge[i]);
    }
    free(lattice->edge);
    free(lattice->adjacent_edge);
    free(lattice);
    return NULL;
  }
  for (i = 0; i < lattice->n_sites; i++) {
    lattice->incident_edge[i] = (int *)malloc(2 * lattice->dimension * sizeof(int));
    if (!lattice->incident_edge[i]) {
      fprintf(stderr, "Memory allocation failed for incident-edge array row %d.\n", i);
      for (j = 0; j < lattice->n_sites; j++) {
      	free(lattice->site[j]);
      	free(lattice->nearest_neighbor[j]);
      }
      free(lattice->site);
      free(lattice->nearest_neighbor);
      for (i = 0; i < lattice->n_edges; i++) {
      	free(lattice->edge[i]);
      	free(lattice->adjacent_edge[i]);
      }
      free(lattice->edge);
      free(lattice->adjacent_edge);
      for (j = 0; j < i; j++) free(lattice->incident_edge[j]);
      free(lattice->incident_edge);
      free(lattice);
      return NULL;
    }
  }

  // READ THE FILE
  // Read sites
  if (fscanf(file, "\n#lattice sites") == EOF) {
    printf("Fatal error while reading input file.\n");
    exit(1);
  }
  for (i = 0; i < lattice->n_sites; i++) {
    if (fscanf(file, "%d", &tmp) == EOF) {
      printf("Fatal error while reading input file.\n");
      exit(1);
    }
    for (j = 0; j < lattice->dimension; j++) {
      if (fscanf(file, "%d", &lattice->site[i][j]) == EOF) {
        printf("Fatal error while reading input file.\n");
        exit(1);
      }
    }
  }

  // Read edges
  if (fscanf(file, "\n#lattice edges") == EOF) {
    printf("Fatal error while reading input file.\n");
    exit(1);
  }
  for (i = 0; i < lattice->n_edges; i++) {
    if (fscanf(file, "%d %d %d", &tmp, &lattice->edge[i][0], &lattice->edge[i][1]) == EOF) {
      printf("Fatal error while reading input file.\n");
      exit(1);
    }
  }
  
  // Read nearest-neighbors
  if (fscanf(file, "\n#nearest-neighbors") == EOF) {
    printf("Fatal error while reading input file.\n");
    exit(1);
  }
  for (i = 0; i < lattice->n_sites; i++) {
    if (fscanf(file, "%d", &tmp) == EOF) {
      printf("Fatal error while reading input file.\n");
      exit(1);
    }
    for (j = 0; j < 2 * lattice->dimension; j++) {
      if (fscanf(file, "%d", &lattice->nearest_neighbor[i][j]) == EOF) {
        printf("Fatal error while reading input file.\n");
        exit(1);
      }
    }
  }

  // Read adjacent edges
  if (fscanf(file, "\n#adjacent edges") == EOF) {
    printf("Fatal error while reading input file.\n");
    exit(1);
  }
  for (i = 0; i < lattice->n_edges; i++) {
    if (fscanf(file, "%d", &tmp) == EOF) {
      printf("Fatal error while reading input file.\n");
      exit(1);
    }
    for (j = 0; j < 2 * (2 * lattice->dimension - 1); j++) {
      if (fscanf(file, "%d", &lattice->adjacent_edge[i][j]) == EOF) {
        printf("Fatal error while reading input file.\n");
        exit(1);
      }
    }
  }

  // Read incident edges
  if (fscanf(file, "\n#incident edges") == EOF) {
    printf("Fatal error while reading input file.\n");
    exit(1);
  }
  for (i = 0; i < lattice->n_sites; i++) {
    if (fscanf(file, "%d", &tmp) == EOF) {
      printf("Fatal error while reading input file.\n");
      exit(1);
    }
    for (j = 0; j < 2 * lattice->dimension; j++) {
      if (fscanf(file, "%d", &lattice->incident_edge[i][j]) == EOF) {
        printf("Fatal error while reading input file.\n");
        exit(1);
      }
    }
  }

  fclose(file);

  return lattice;
}


// HYPER_CUBIC_LATTICE DESTROYER
void free_lattice(HC_Lattice *lattice) {
  if (lattice) {
    int i;
    for (i = 0; i < lattice->n_sites; i++) {
      free(lattice->site[i]);
      free(lattice->nearest_neighbor[i]);
      free(lattice->incident_edge[i]);
    }
    free(lattice->site);
    free(lattice->nearest_neighbor);
    free(lattice->incident_edge);
    
    for (i = 0; i < lattice->n_edges; i++) {
      free(lattice->edge[i]);
      free(lattice->adjacent_edge[i]);
    }
    free(lattice->edge);
    free(lattice->adjacent_edge);
    
    free(lattice);
  }
}


//////////////////////////////////////////////////////////////////////////////
// FUNTION TO EXTRACT THE LINEAR CHAIN
Chain *extract_linear_chain(HC_Lattice *lattice, char *edge_string, int endpoint) {
  int i, j, curr_bond, counter, start_site, next_site;
  
  start_site = endpoint;
  for (i = 0; i < 2 * lattice->dimension; i++) {
    curr_bond = lattice->incident_edge[start_site][i];
    if (edge_string[curr_bond] == '1') break;
  }
  
  if (start_site == lattice->edge[curr_bond][0]) next_site = lattice->edge[curr_bond][1];
  else next_site = lattice->edge[curr_bond][0];

  // Initialize chain array
	int length = 2;
  int *chain_array = (int *)malloc(length * sizeof(int));

  chain_array[0] = start_site;
	chain_array[1] = next_site;

	double mean_sq_displ = 0.0;
	double *Rcm = (double *)malloc(lattice->dimension * sizeof(double));

	int *cell = (int *)malloc(lattice->dimension * sizeof(int));
	for (i = 0; i < lattice->dimension; i++) {
		if ((lattice->site[start_site][i] == 0) && (lattice->site[next_site][i] == lattice->size - 1)) cell[i] = -1;
		else if ((lattice->site[start_site][i] == lattice->size - 1) && (lattice->site[next_site][i] == 0)) cell[i] = 1;
		else cell[i] = 0;
		Rcm[i] = (lattice->site[start_site][i] + lattice->site[next_site][i] + lattice->size * cell[i]) / 2.0;
		mean_sq_displ += (lattice->site[start_site][i] + lattice->size * cell[i]) * (lattice->site[start_site][i] + lattice->size * cell[i]) + (lattice->site[next_site][i] + lattice->size * cell[i]) * (lattice->site[next_site][i] + lattice->size * cell[i]); 
	}
	mean_sq_displ /= 2.0;

  // Initialize coordinates matrix
  int **chain_coordinate = (int **)malloc(length * sizeof(int *));
  for (i = 0; i < length; i++) {
    chain_coordinate[i] = (int *)malloc(lattice->dimension * sizeof(int));
    if (!chain_coordinate[i]) {
      fprintf(stderr, "Memory allocation failed for chain_coordinate of %d-th chain monomer.\n", i);
      for (j = 0; j < i; j++) free(chain_coordinate[i]);
      free(chain_array);
    }
  }

  for (i = 0; i < lattice->dimension; i++) {
    chain_coordinate[0][i] = lattice->site[start_site][i];
    chain_coordinate[1][i] = lattice->site[next_site][i] + lattice->size * cell[i];
  }

	start_site = next_site;

	// Traverse the chain 
	counter = 2;
  while (next_site != -1) {
		next_site = find_next_site(lattice, edge_string, start_site, &curr_bond);
		if (next_site == -1) break;

		mean_sq_displ *= counter;
		for (i = 0; i < lattice->dimension; i++) {
			if ((lattice->site[start_site][i] == 0) && (lattice->site[next_site][i] == lattice->size - 1)) cell[i]--;
			if ((lattice->site[start_site][i] == lattice->size - 1) && (lattice->site[next_site][i] == 0)) cell[i]++;
			Rcm[i] = (counter * Rcm[i] + lattice->site[next_site][i] + lattice->size * cell[i]) / (counter + 1);
			mean_sq_displ += (lattice->site[next_site][i] + lattice->size * cell[i]) * (lattice->site[next_site][i] + lattice->size * cell[i]);
		}
		counter++;
		mean_sq_displ /= counter;

		start_site = next_site;
		append_element_to_int_array(&chain_array, &length, next_site);

    // Reallocate coordinates
    chain_coordinate = (int **)realloc(chain_coordinate, length * sizeof(int *));
    chain_coordinate[length - 1] = (int *)malloc(lattice->dimension * sizeof(int));
		for (i = 0; i < lattice->dimension; i++) chain_coordinate[length - 1][i] = lattice->site[next_site][i] + lattice->size * cell[i];
	}

	length--;

	// Compute Rg and Ree
	double dist_sq = 0.0;
	for (i = 0; i < lattice->dimension; i++) {
		mean_sq_displ -= Rcm[i] * Rcm[i];
		dist_sq += (lattice->site[chain_array[length]][i] + lattice->size * cell[i] - lattice->site[chain_array[0]][i]) * (lattice->site[chain_array[length]][i] + lattice->size * cell[i] - lattice->site[chain_array[0]][i]); 
	}
	double Rg = sqrt(mean_sq_displ);
	double Ree = sqrt(dist_sq);

  // Create the chain
  Chain *chain = (Chain *)malloc(sizeof(Chain));
  chain->linear = true;
  chain->length = length;

  chain->chain_array = (int *)malloc((length + 1) * sizeof(int));
	for (i = 0; i <= length; i++) chain->chain_array[i] = chain_array[i]; 

  chain->chain_coordinate = (int **)malloc((length + 1) * sizeof(int *));
  for (i = 0; i <= length; i++) {
    chain->chain_coordinate[i] = (int *)malloc(lattice->dimension * sizeof(int));
    for (j = 0; j < lattice->dimension; j++) chain->chain_coordinate[i][j] = chain_coordinate[i][j];
  }

  chain->Rcm = (double *)malloc(lattice->dimension * sizeof(double));
  for (i = 0; i < lattice->dimension; i++) chain->Rcm[i] = Rcm[i];

	chain->Rg = Rg;
	chain->Ree = Ree;

  // Free memory
  free(chain_array);
  free(Rcm);
  free(cell);
  for (i = 0; i <= length; i++) free(chain_coordinate[i]);
  free(chain_coordinate);

  return chain;
}


//////////////////////////////////////////////////////////////////////////////
// FUNCTION TO EXTRACT A RING
Chain *extract_ring(HC_Lattice *lattice, char *edge_string) {
  int i, counter, curr_bond, start_site, next_site;
	double mean_sq_displ = 0.0;
	double *Rcm = (double *)malloc(lattice->dimension * sizeof(double));

  // Find the start of the ring
  i = 0;
	while ((i < lattice->n_edges) && (edge_string[i] == '0')) i++;
	if (i == lattice->n_edges) return NULL; // no active edge found

	// Initialize ring array
	start_site = lattice->edge[i][0];
	next_site = lattice->edge[i][1];
	curr_bond = i;

	int length = 2;
	int *ring_array = (int *)malloc(length * sizeof(int));
	ring_array[0] = start_site;
	ring_array[1] = next_site;

	edge_string[i] = '0'; // remove the bond from the configuration

	int *cell = (int *)malloc(lattice->dimension * sizeof(int));
	for (i = 0; i < lattice->dimension; i++) {
		if ((lattice->site[start_site][i] == 0) && (lattice->site[next_site][i] == lattice->size - 1)) cell[i] = -1;
		else if ((lattice->site[start_site][i] == lattice->size - 1) && (lattice->site[next_site][i] == 0)) cell[i] = 1;
		else cell[i] = 0;
	}

	for (i = 0; i < lattice->dimension; i++) {
		Rcm[i] = (lattice->site[ring_array[0]][i] + lattice->site[ring_array[1]][i] + lattice->size * cell[i]) / 2.0;
		mean_sq_displ += lattice->site[ring_array[0]][i] * lattice->site[ring_array[0]][i] + (lattice->site[ring_array[1]][i] + lattice->size * cell[i]) * (lattice->site[ring_array[1]][i] + lattice->size * cell[i]); 
	}
	mean_sq_displ /= 2.0;

  int **ring_coordinate = (int **)malloc(length * sizeof(int *));
  for (i = 0; i < length; i++) {
    ring_coordinate[i] = (int *)malloc(lattice->dimension * sizeof(int));
    if (!ring_coordinate[i]) {
      fprintf(stderr, "Memory allocation failed for ring_coordinate of %d-th monomer.\n", i);
      int j;
      for (j = 0; j < i; j++) free(ring_coordinate[j]);
      free(ring_array);
    }
  }

  for (i = 0; i < lattice->dimension; i++) {
    ring_coordinate[0][i] = lattice->site[start_site][i];
    ring_coordinate[1][i] = lattice->site[next_site][i] + lattice->size * cell[i];
  }

	for (i = 0; i < 2 * (2 * lattice->dimension - 1); i++) {
		counter = lattice->adjacent_edge[curr_bond][i];
		if (edge_string[counter] == '1') {
			if (lattice->edge[counter][0] == next_site)	{
				next_site = lattice->edge[counter][1];
				break;
			}	else if (lattice->edge[counter][1] == next_site) {
				next_site = lattice->edge[counter][0];
				break;
			}
		}
	}
	append_element_to_int_array(&ring_array, &length, next_site);

	mean_sq_displ *= 2.0;
	for (i = 0; i < lattice->dimension; i++) {
		if ((lattice->site[ring_array[1]][i] == 0) && (lattice->site[next_site][i] == lattice->size - 1)) cell[i]--;
		if ((lattice->site[ring_array[1]][i] == lattice->size - 1) && (lattice->site[next_site][i] == 0)) cell[i]++;
		Rcm[i] = (2 * Rcm[i] + lattice->site[next_site][i] + lattice->size * cell[i]) / 3.0;
		mean_sq_displ += (lattice->site[next_site][i] + lattice->size * cell[i]) * (lattice->site[next_site][i] + lattice->size * cell[i]);
	}
	mean_sq_displ /= 3.0;

  // Reallocate coordinates
  int **new_coordinate = (int **)realloc(ring_coordinate, length * sizeof(int *));
  ring_coordinate = new_coordinate;
  ring_coordinate[length - 1] = (int *)malloc(lattice->dimension * sizeof(int));
  for (i = 0; i < lattice->dimension; i++) ring_coordinate[length - 1][i] = lattice->site[next_site][i] + lattice->size * cell[i];

	curr_bond = counter;
	start_site = next_site;
	edge_string[curr_bond] = '0';

	// Traverse the ring 
	counter = 3;
  while (next_site != ring_array[0]) {
		next_site = find_next_site(lattice, edge_string, start_site, &curr_bond);
		if (next_site == ring_array[0]) {
			edge_string[curr_bond] = '0';
			for (i = 0; i < lattice->dimension; i++) {
				if ((lattice->site[start_site][i] == 0) && (lattice->site[next_site][i] == lattice->size - 1)) cell[i]--;
				if ((lattice->site[start_site][i] == lattice->size - 1) && (lattice->site[next_site][i] == 0)) cell[i]++;
			}
			break;
		}

		mean_sq_displ *= counter;
		for (i = 0; i < lattice->dimension; i++) {
			if ((lattice->site[start_site][i] == 0) && (lattice->site[next_site][i] == lattice->size - 1)) cell[i]--;
			if ((lattice->site[start_site][i] == lattice->size - 1) && (lattice->site[next_site][i] == 0)) cell[i]++;
			Rcm[i] = (counter * Rcm[i] + lattice->site[next_site][i] + lattice->size * cell[i]) / (counter + 1);
			mean_sq_displ += (lattice->site[next_site][i] + lattice->size * cell[i]) * (lattice->site[next_site][i] + lattice->size * cell[i]);
		}
		counter++;
		mean_sq_displ /= counter;

		start_site = next_site;
		append_element_to_int_array(&ring_array, &length, next_site);

    // Reallocate coordinates
    int **new_coordinate = (int **)realloc(ring_coordinate, length * sizeof(int *));
    ring_coordinate = new_coordinate;
    ring_coordinate[length - 1] = (int *)malloc(lattice->dimension * sizeof(int));
		for (i = 0; i < lattice->dimension; i++) ring_coordinate[length - 1][i] = lattice->site[next_site][i] + lattice->size * cell[i];
	}

	// Compute Rg
	double Rg;
	counter = 0;
	for (i = 0; i < lattice->dimension; i++) {
		mean_sq_displ -= Rcm[i] * Rcm[i];
		counter += abs(cell[i]);
	}
	if (counter == 0) Rg = sqrt(mean_sq_displ);
	else Rg = -1.0; // set Rg = -1 if the ring is unphysical (i.e. it closes only on the thorus and does not in the periodic lattice)

  // Create the chain
  Chain *ring = (Chain *)malloc(sizeof(Chain));
  ring->linear = false;
  ring->length = length;

  ring->chain_array = (int *)malloc(length * sizeof(int));
	for (i = 0; i < length; i++) ring->chain_array[i] = ring_array[i]; 

  ring->chain_coordinate = (int **)malloc(length * sizeof(int *));
	int j;
  for (i = 0; i < length; i++) {
    ring->chain_coordinate[i] = (int *)malloc(lattice->dimension * sizeof(int));
    for (j = 0; j < lattice->dimension; j++) ring->chain_coordinate[i][j] = ring_coordinate[i][j];
  }

  ring->Rcm = (double *)malloc(lattice->dimension * sizeof(double));
  for (i = 0; i < lattice->dimension; i++) ring->Rcm[i] = Rcm[i];

	ring->Rg = Rg;
	ring->Ree = -1;

  // Free memory
  free(ring_array);
  free(Rcm);
  free(cell);
  for (i = 0; i < length; i++) free(ring_coordinate[i]);
  free(ring_coordinate);

  return ring;
}


//////////////////////////////////////////////////////////////////////////////
// FUNCTION TO FIND NEXT SITE OF A CHAIN (return -1 if no site has been found)
int find_next_site(HC_Lattice *lattice, char *edge_string, int start_site, int *curr_bond) {
	edge_string[(*curr_bond)] = '0';
	int i;
	bool next_bond_found = false;
	for (i = 0; i < 2 * (2 * lattice->dimension - 1); i++) {
		if (edge_string[lattice->adjacent_edge[(*curr_bond)][i]] == '1') {
			next_bond_found = true;
			break;
		}
	}
	(*curr_bond) = lattice->adjacent_edge[(*curr_bond)][i];
	
	if (!next_bond_found) return -1;
	if (lattice->edge[(*curr_bond)][0] == start_site) return lattice->edge[(*curr_bond)][1];
	else return lattice->edge[(*curr_bond)][0];
}


//////////////////////////////////////////////////////////////////////////////
// FUNCTION TO APPEND AN ELEMENT TO AN INT ARRAY
void append_element_to_int_array(int **array, int *size, int element) {
  // Reallocate memory to increase the size of the array
  *array = (int *)realloc(*array, (*size + 1) * sizeof(int));
  if (!(*array)) {
    printf("Memory allocation failed!\n");
    exit(1);
  }
  (*array)[*size] = element;
  (*size)++;
}


//////////////////////////////////////////////////////////////////////////////
// CHAIN DESTROYER
void free_chain(Chain *chain) {
  if (chain) {
	  if (chain->chain_array) free(chain->chain_array);
	  if (chain->Rcm) free(chain->Rcm);
	  if (chain->chain_coordinate) {
	    int i;
	    int chain_array_size = chain->length;
	    if (chain->linear) chain_array_size++;
	    for (i = 0; i < chain_array_size; i++) if (chain->chain_coordinate[i]) free(chain->chain_coordinate[i]);
	    free(chain->chain_coordinate);
	  }
  	free(chain);
	}
}


//////////////////////////////////////////////////////////////////////////////
// FUNCTION TO EXTRACT THE CHAINS FROM THE BOND STRING
Configuration *construct_configuration_chains(HC_Lattice *lattice, char *edge_string, int *end, int n_linear_chains) {
  Configuration *config = (Configuration *)malloc(sizeof(Configuration));
  if (!config) {
    perror("Failed to allocate memory for Configuration");
    return NULL;
  }

  config->n_rings = 0;
	config->n_linear_chains = n_linear_chains;
  config->ring = NULL;

	int i, j, k;
  config->linear_chain = (Chain **)malloc(n_linear_chains * sizeof(Chain *));
	for (i = 0; i < n_linear_chains; i++) config->linear_chain[i] = NULL;

  // Find the linear chains
	int *aux_end = (int *)malloc(2 * n_linear_chains * sizeof(int));
	for (i = 0; i < 2 * n_linear_chains; i++) aux_end[i] = end[i];

	int aux_endpoint;
	for (i = 0; i < n_linear_chains; i++) {
		j = 0;
		while (aux_end[j] == -2) j++;
		config->linear_chain[i] = extract_linear_chain(lattice, edge_string, aux_end[j]);
		aux_end[j] = -2;
		aux_endpoint = config->linear_chain[i]->chain_array[config->linear_chain[i]->length];
		for (k = 0; k < 2 * n_linear_chains; k++) if (aux_end[k] == aux_endpoint) aux_end[k] = -2;
	}
	free(aux_end);

	// Choose linear_chain[0] having center-of-mass in lattice cell 0
  int replica;
	for (i = 0; i < lattice->dimension; i++) {
		replica = (int)round((double)(config->linear_chain[0]->Rcm[i]) / lattice->size);
		config->linear_chain[0]->Rcm[i] -= lattice->size * replica;
		for (j = 0; j <= config->linear_chain[0]->length; j++) config->linear_chain[0]->chain_coordinate[j][i] -= lattice->size * replica;
	}

	// Choose other linear chains replicas as close to linear_chain[0] as possible
	if (n_linear_chains > 1) {
		for (k = 1; k < n_linear_chains; k++) {
			for (i = 0; i < lattice->dimension; i++) {
				replica = (int)round((double)(config->linear_chain[k]->Rcm[i] - config->linear_chain[0]->Rcm[i]) / lattice->size);
				config->linear_chain[k]->Rcm[i] -= lattice->size * replica;
				for (j = 0; j <= config->linear_chain[k]->length; j++) config->linear_chain[k]->chain_coordinate[j][i] -= lattice->size * replica;
			}
		}
	}

  // Extract all rings until there are no more
  Chain *ring;
  while (ring = extract_ring(lattice, edge_string)) {
 		config->ring = (Chain **)realloc(config->ring, (config->n_rings + 1) * sizeof(Chain *));
  	if (!config->ring) {
      perror("Failed to allocate memory for rings");
      free_configuration(config);
      return NULL;
    }
    config->ring[config->n_rings] = ring;
    config->n_rings++;
  }

  return config;
}


// CONFIGURATION DESTROYER
void free_configuration(Configuration *config) {
	if (config) {
		int i;
		if (config->linear_chain) {
	    for (i = 0; i < config->n_linear_chains; i++) if (config->linear_chain[i]) free_chain(config->linear_chain[i]);
  	  free(config->linear_chain);
		}
		if (config->ring) {
			for (i = 0; i < config->n_rings; i++) if (config->ring[i]) free_chain(config->ring[i]);
			free(config->ring);
		}
		free(config);
	}
}
