///////////////////////////////////////////////////////////
// efornasa, February 2026                               //
// Implementation of the functions defined in the header //
// file "functions.h" for SAMC sampling.                 //
///////////////////////////////////////////////////////////
// FUNCTIONS:                                            //
//                                                       //
// construct_paramenters................................ //
// construct_lattice.................................... //
// construct_lattice_from_input_file.................... //
// free_lattice......................................... //
// initialize_configuration............................. //
// initialize_configuration_from_input_file............. //
// free_configuration................................... //
// metropolis_move...................................... //
//                                                       //
// EVOLUTION............................................ //
//                                                       //
// RAN2 RANDOM GENERATOR................................ //
/////////////////////////////////////////////////////////// 

// functions.c

#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdbool.h>
#include <time.h>

#include "functions.h"
#include "ran2.c" // Random number generator from "Numerical Recipes" by Press et al. Users must obtain their own licensed version


//////////////////////////////////////////////////////////////////////////////
// PARAMETERS CONSTRUCTOR
Parameters *construct_parameters(const char *input_file) {
  FILE *file = fopen(input_file, "r");
  if (!file) {
    printf("Error: Could not open file %s\n", input_file);
    return NULL;
  }
  
  Parameters *param = (Parameters *)malloc(sizeof(Parameters));
  if (!param) {
    printf("Error: Could not allocate memory for parameters.\n");
    fclose(file);
    return NULL;
  }
  
  char *key = (char *)malloc(64 * sizeof(char));
  int value;
  
  while (fscanf(file, "%s %d", key, &value) == 2) {
    if (strcmp(key, "dimension") == 0) param->dimension = value;
    else if (strcmp(key, "size") == 0) param->size = value;
    else if (strcmp(key, "n_linear_chains") == 0) param->n_linear_chains = value;
    else if (strcmp(key, "n_dumps") == 0) param->n_dumps = value;
    else if (strcmp(key, "n_steps_per_dump") == 0) param->n_steps_per_dump = value;
  }
  fclose(file);
	free(key);

  return param;
}


//////////////////////////////////////////////////////////////////////////////
// LATTICE CONSTRUCTOR
HC_Lattice *construct_lattice(int dimension, int size) {
  // Allocate memory for the HC_Lattice structure
  HC_Lattice *lattice = (HC_Lattice *)malloc(sizeof(HC_Lattice));
  if (!lattice) {
    fprintf(stderr, "Memory allocation failed for HC_Lattice.\n");
    return NULL;
  }
  
  int i, j;

  // Initialize the dimension and size
  lattice->dimension = dimension;
  lattice->size = size;

  // Set the number of sites and edges
  lattice->n_sites = pow(size, dimension);

  lattice->n_edges = dimension * lattice->n_sites;
 
  // MEMORY ALLOCATION  
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
      // Free previously allocated memory before returning
      for (j = 0; j < i; j++) free(lattice->site[j]);
      free(lattice->site);
      free(lattice);
      return NULL;
    }
  }
  printf("\nLATTICE_CONSTRUCTOR\nMemory allocated for site array\n");
  
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
    if (!lattice->nearest_neighbor[i]) {
      fprintf(stderr, "Memory allocation failed for nearest-neighbor array row %d.\n", i);
      // Free previously allocated memory before returning
      for (j = 0; j < lattice->n_sites; j++) free(lattice->site[j]);
      free(lattice->site);
      for (j = 0; j < i; j++) free(lattice->nearest_neighbor[j]);
      free(lattice->nearest_neighbor);
      free(lattice);
      return NULL;
    }
  }
  printf("Memory allocated for nearest-neighbor array\n");
  
  // Allocate memory for edge array
  lattice->edge = (int **)malloc(lattice->n_edges * sizeof(int *));
  if (!lattice->edge) {
    fprintf(stderr, "Memory allocation failed for edge array.\n");
    // Free previously allocated memory before returning
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
      // Free previously allocated memory before returning
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
  printf("Memory allocated for edge array\n");
  
  // Allocate memory for adjacent-edge array
  lattice->adjacent_edge = (int **)malloc(lattice->n_edges * sizeof(int *));
  if (!lattice->adjacent_edge) {
    fprintf(stderr, "Memory allocation failed for adjacent-edge array.\n");
    // Free previously allocated memory before returning
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
      fprintf(stderr, "Memory allocation failed for adjacent-edge array row %d.\n", i);
      // Free previously allocated memory before returning
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
  printf("Memory allocated for adjacent-edge array\n");
  
  // Allocate memory for incident-edge array
  lattice->incident_edge = (int **)malloc(lattice->n_sites * sizeof(int *));
  if (!lattice->incident_edge) {
    fprintf(stderr, "Memory allocation failed for incident-edge array.\n");
    // Free previously allocated memory before returning
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
      // Free previously allocated memory before returning
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
  printf("Memory allocated for incident-edge array\n\n");
  
  // FILL ARRAYS: set to -1 elements that do not exist
  // Populate site array
  int index, aux;
  for (i = 0; i < lattice->n_sites; i++) {
    index = i;
    aux = lattice->n_sites;
    for (j = dimension - 1; j >= 1; j--) {
      aux /= size;
      lattice->site[i][j] = index / aux;
      index -= lattice->site[i][j] * aux;
    }
    lattice->site[i][0] = index;
  }
  printf("site array filled\n");
 
  // Populate nearest-neighbor array
  for (i = 0; i < lattice->n_sites; i++) {
    aux = 1;
    for (j = 0; j < dimension; j++) {
      lattice->nearest_neighbor[i][2 * j] = i - aux;
      lattice->nearest_neighbor[i][2 * j + 1] = i + aux;
      aux *= size;
      if (lattice->site[i][j] == 0) lattice->nearest_neighbor[i][2 * j] += aux;
      else if (lattice->site[i][j] == (lattice->size - 1)) lattice->nearest_neighbor[i][2 * j + 1] -= aux;
    }
  }	
  printf("nearest-neighbor array filled\n");
  
  // Populate edge and incident-edge array
  aux = 0;
  int *counter = (int *)malloc(lattice->n_sites * sizeof(int));
  for (i = 0; i < lattice->n_sites; i++) counter[i] = 0;
  for (i = 0; i < lattice->n_sites; i++) {
    for (j = 0; j < 2 * dimension; j++) {
      index = lattice->nearest_neighbor[i][j];
      if (index == -1) {
        lattice->incident_edge[i][counter[i]] = -1;
        counter[i]++;
      }
      else if (index > i) {//printf("NN of site %d: %d\n", i, index);
      	lattice->edge[aux][0] = i;
      	lattice->edge[aux][1] = index;
        lattice->incident_edge[i][counter[i]] = aux;
        lattice->incident_edge[index][counter[index]] = aux;
        counter[i]++;
        counter[index]++;
      	aux++;
      }
    }
  }
  free(counter);
  printf("edge and incident-edge arrays filled\n");

  // Populate adjacent-edge array
  for (i = 0; i < lattice->n_edges; i++) {
    aux = 0;
    index = lattice->edge[i][0];//printf("edge %d, site %d\n", i, index);
    for (j = 0; j < 2 * dimension; j++) {
      if (lattice->incident_edge[index][j] == -1) {
        lattice->adjacent_edge[i][aux] = -1;
        aux++;
      }
      else if (lattice->incident_edge[index][j] != i) {
        lattice->adjacent_edge[i][aux] = lattice->incident_edge[index][j];
        aux++;
      }
    }
    index = lattice->edge[i][1];
    for (j = 0; j < 2 * dimension; j++) {
      if (lattice->incident_edge[index][j] == -1) {
        lattice->adjacent_edge[i][aux] = -1;
        aux++;
      }
      else if (lattice->incident_edge[index][j] != i) {
        lattice->adjacent_edge[i][aux] = lattice->incident_edge[index][j];
        aux++;
      }
    }
  } 
  printf("adjacent-edge array filled\n\n");

 
  return lattice;
}


//////////////////////////////////////////////////////////////////////////////
// HYPER_CUBIC_LATTICE CONSTRUCTOR FROM INPUT FILE
HC_Lattice *construct_lattice_from_input_file(const char *input_lattice_filename) {
  FILE *file = fopen(input_lattice_filename, "r");
  if (!file) {
    printf("Error: Could not open file %s\n", input_lattice_filename);
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
  if(fscanf(file, "#size %d\n", &lattice->size) == EOF) {
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
      // Free previously allocated memory before returning
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
      // Free previously allocated memory before returning
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
    // Free previously allocated memory before returning
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
      // Free previously allocated memory before returning
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
    // Free previously allocated memory before returning
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
      // Free previously allocated memory before returning
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
    // Free previously allocated memory before returning
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
      // Free previously allocated memory before returning
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

  for (i = 0; i < lattice->n_edges; i++){
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

    for (j = 0; j < 2 * lattice->dimension; j++){
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
    for (j = 0; j < 2 * (2 * lattice->dimension - 1); j++){
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
    for (j = 0; j < 2 * lattice->dimension; j++){
      if (fscanf(file, "%d", &lattice->incident_edge[i][j]) == EOF) {
        printf("Fatal error while reading input file.\n");
        exit(1);
      }
    }
  }

  fclose(file);

  return lattice;
}


//////////////////////////////////////////////////////////////////////////////
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
// INITIALIZE CONFIGURATION
Configuration *initialize_configuration(HC_Lattice *lattice, int n_linear_chains) {
  Configuration *config = (Configuration *)malloc(sizeof(Configuration));
  if (!config) {
    perror("Failed to allocate memory for the configuration\n");
    return NULL;
  }
  
  config->n_spins = lattice->n_edges;
	config->n_linear_chains = n_linear_chains;

	config->end = (int *)malloc(2 * n_linear_chains * sizeof(int));
	if(!config->end) {
    perror("Failed to allocate memory for the ends\n");
    free(config);
    return NULL;
  }
  
  config->spin = (char *)malloc(config->n_spins * sizeof(char));
  if (!config->spin) {
    perror("Failed to allocate memory for the spins\n");
		free(config->end);
    free(config);
    return NULL;
  }

  // Set to '1' all horizontal edges
  int i, j, k, e;
  for (i = 0; i < lattice->n_edges; i++) {
    if (lattice->site[lattice->edge[i][0]][0] == (lattice->site[lattice->edge[i][1]][0] + 1) % lattice->size) config->spin[i] = '1';
    else if (lattice->site[lattice->edge[i][1]][0] == (lattice->site[lattice->edge[i][0]][0] + 1) % lattice->size) config->spin[i] = '1';
    else config->spin[i] = '0';
  }

  // Chose deterministically the endpoints and the n_linear_chains horizontal edges to be set to '0'
	e = 0;
	for (i = 0; i < n_linear_chains; i++) {
		if (e % lattice->size == 0) e += lattice->size; 
		if ((e + 1) % lattice->size == 0) e += lattice->size + 1; 
		config->end[i] = e;
		config->end[i + n_linear_chains] = e + 1;
		for (j = 0; j < 2 * lattice->dimension; j++) {
			k = lattice->incident_edge[e][j];
			if (config->spin[k] == '1') if (((lattice->edge[k][0] == e) && (lattice->edge[k][1] == e + 1)) || ((lattice->edge[k][0] == e + 1) && (lattice->edge[k][1] == e))) config->spin[k] = '0';
		}
		e += 2;
	}

  return config;
}


//////////////////////////////////////////////////////////////////////////////
// INITIALIZE CONFIGURATION FROM INPUT FILE
Configuration *initialize_configuration_from_input_file(const char *input_file, int n_spins, int n_linear_chains) {
  FILE *file = fopen(input_file, "r");
  if (!file) {
    printf("Error: Could not open file %s\n", input_file);
    return NULL;
  }

  Configuration *config = malloc(sizeof(Configuration));
  if (!config) {
    perror("Failed to allocate memory for the configuration\n");
    fclose(file);
    return NULL;
  }

	config->n_spins = n_spins;
	config->n_linear_chains = n_linear_chains;

  config->spin = (char *)malloc(n_spins * sizeof(char));
  if (!config->spin) {
    perror("Failed to allocate memory for spins");
		free(config);
    fclose(file);
    return NULL;
  }

  config->end = (int *)malloc(2 * n_linear_chains * sizeof(int));
  if (!config->spin) {
    perror("Failed to allocate memory for ends");
		free(config->spin);
		free(config);
    fclose(file);
    return NULL;
  }

  if (fscanf(file, "%s", config->spin) == EOF) {
    printf("Fatal error while reading input file.\n");
    exit(1);
  }

  int i;
  for (i = 0; i < 2 * n_linear_chains; i++) {
    if (fscanf(file, "%d", &config->end[i]) == EOF) {
      printf("Fatal error while reading input file.\n");
      exit(1);
    }
  }
  fclose(file);  

  return config;
}


//////////////////////////////////////////////////////////////////////////////
// CONFIGURATION DESTROYER
void free_configuration(Configuration *config) {
  if (config) {
    if (config->end) free(config->end);
    if (config->spin) free(config->spin);
    free(config);
  }
}


//////////////////////////////////////////////////////////////////////////////
// METROPOLIS MOVE
void metropolis_move(Configuration *config, HC_Lattice *lattice, long *idum) {
  int e, i, j, site_A, site_B, site_C, edge_index, new_bond, old_bond;
  int n_incident_edges = 2 * lattice->dimension;
	int n_ends = 2 * config->n_linear_chains;

	// Select randomly from which end of the open chains has to change (site_A)
	e = rand_int(n_ends, idum);
	site_A = config->end[e];
  
  // Find the end-bond (incident_edge[A][i])
  for (i = 0; i < n_incident_edges; i++) {
    edge_index = lattice->incident_edge[site_A][i];
		if (config->spin[edge_index] == '1') break;
  }

  // Choose which bond activate (incident_edge[A][j])
	j = rand_int(n_incident_edges - 1, idum);
  if (j >= i) j++;
	new_bond = lattice->incident_edge[site_A][j];

	// Find site B that makes the new bond with A
	if (lattice->edge[new_bond][0] == site_A) site_B = lattice->edge[new_bond][1];
	else site_B = lattice->edge[new_bond][0];

	// Reject the move if site_B is one of the ends
  for (i = 0; i < n_ends; i++) if (site_B == config->end[i]) return;

  // Choose which bond inactivate (old_bond)
	j = 0; // counts bonds inciding site_B
	for (i = 0; i < n_incident_edges; i++) {
		edge_index = lattice->incident_edge[site_B][i];
		if (config->spin[edge_index] == '1') {
			old_bond = edge_index;
			if (j == 0) {
				if (ran2(idum) < 0.5) break;
				else j++;
			}
			else break;
		}
	}

	if (site_B == lattice->edge[old_bond][0]) site_C = lattice->edge[old_bond][1];
	else site_C = lattice->edge[old_bond][0];

	// Reject if C is already endpoint
	for (i = 0; i < n_ends; i++) if (site_C == config->end[i]) return;

	// Reset the configuration
  config->spin[new_bond] = '1';
  config->spin[old_bond] = '0';
	config->end[e] = site_C;
}


//////////////////////////////////////////////////////////////////////////////
// EVOLUTION
void EVOLUTION(Configuration *config, HC_Lattice *lattice, int n_dumps, int n_steps_per_dump, long *idum) {
  int i, j;
  
  char *system = (char *)malloc(64 * sizeof(char));
  if (lattice->dimension == 2) snprintf(system, 64, "%dx%d", lattice->size, lattice->size);
  else if (lattice->dimension == 3) snprintf(system, 64, "%dx%dx%d", lattice->size, lattice->size, lattice->size);
  else if (lattice->dimension == 4) snprintf(system, 64, "%dx%dx%dx%d", lattice->size, lattice->size, lattice->size, lattice->size);
  else {
    printf("Error: Unsupported lattice dimension %d\n", lattice->dimension);
    return;
  }
  
  // Trajectory
	clock_t start, now;
	double cpu_time;
	start = clock();

	char *rand_state_filename = (char *)malloc(64 * sizeof(char));
	snprintf(rand_state_filename, 64, "log/log_rand_state.txt");
	FILE *rand_state_file = fopen(rand_state_filename, "w");

  for (i = 0; i < n_dumps; i++) {
    for (j = 0; j < n_steps_per_dump; j++) metropolis_move(config, lattice, idum); // Evolve n_steps_per_dump times

		// Save the random number state (idum)
		fprintf(rand_state_file, "%ld\n", *idum);

  	char *filename = (char *)malloc(256 * sizeof(char));
		snprintf(filename, 256, "configurations/config_%s_%06d.dat", system, i);
    FILE *file = fopen(filename, "w");
    if (!file) {
      printf("Error: Could not open file %s for writing.\n", filename);
      return;
    }

    for (j = 0; j < config->n_spins; j++) fprintf(file, "%c", config->spin[j]);
    for (j = 0; j < 2 * config->n_linear_chains; j++) fprintf(file, " %d", config->end[j]);

		free(filename);
  	fclose(file);
    
    // Print on terminal
    if ((i + 1) % (n_dumps / 10) == 0) {
			now = clock();
			cpu_time = ((double)(now - start)) / CLOCKS_PER_SEC;
			printf("%d configurations generated in %g seconds\n", n_dumps / 10, cpu_time);
			start = now;
		}
  }
  printf("End of trajectory. %d configurations have been generated\n", n_dumps);
  
	fclose(rand_state_file);
	free(rand_state_filename);
	free(system);
}


int rand_int(int max_value, long *idum) { return max_value * ran2(idum); }

