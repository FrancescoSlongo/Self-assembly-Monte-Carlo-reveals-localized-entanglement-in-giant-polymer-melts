///////////////////////////////////////////////////////////
// efornasa, March 2026                                  //
// SAMC functions implementation.                        //
///////////////////////////////////////////////////////////
// FUNCTIONS:                                            //
//                                                       //
// generate_paramenters.............................. 45 //
//                                                       //
// get_coordinates_from_site_index................... 78 //
// get_site_index_from_coordinates................... 93 //
// generate_nearest_neighbors_list.................. 108 //
//                                                       //
// build_snake_recursive............................ 138 //
// initialize_configuration......................... 160 //
// copy_configuration_from_config_file.............. 234 //
// free_configuration............................... 279 //
//                                                       //
// SAMC_MOVE........................................ 290 //
// samc_trajectory.................................. 332 //
//                                                       //
// rand_int......................................... 372 //
//                                                       //
// alloc_2d_ull..................................... 377 //
// free_2d_ull...................................... 395 //
//                                                       //
// print_unwrapped_configuration.................... 404 //
// periodic_step.................................... 535 //
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
#include "ran2.h"


//////////////////////////////////////////////////////////////////////////////
// PARAMETERS CONSTRUCTOR
Parameters *generate_parameters(const char *input_file) {
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
  size_t value;
  
  while (fscanf(file, "%s %zu", key, &value) == 2) {
    if (strcmp(key, "dimension") == 0) param->dimension = value;
    else if (strcmp(key, "lattice_side") == 0) param->lattice_side = value;
    else if (strcmp(key, "n_linear_chains") == 0) param->n_linear_chains = value;
    else if (strcmp(key, "n_dumps") == 0) param->n_dumps = value;
    else if (strcmp(key, "n_steps_per_dump") == 0) param->n_steps_per_dump = value;
  }
  fclose(file);
	free(key);

  return param;
}


//////////////////////////////////////////////////////////////////////////////
// GET COORDINATES FROM SITE INDEX
size_t *get_coordinates_from_site_index(size_t dimension, size_t lattice_side, size_t site_index) {
	size_t *x = (size_t *)malloc(dimension * sizeof(size_t));
	size_t i;
	size_t aux = site_index;
	for (i = 0; i < dimension; i++) {
		x[i] = aux % lattice_side;
		aux /= lattice_side;
	}

	return x;
}


//////////////////////////////////////////////////////////////////////////////
// GET COORDINATES FROM SITE INDEX
size_t get_site_index_from_coordinates(size_t dimension, size_t lattice_side, size_t *x) {
	size_t i;
	size_t site_index = 0;
	size_t aux = 1;
	for (i = 0; i < dimension; i++) {
		site_index += aux * x[i];
		aux *= lattice_side;
	}

	return site_index;
}


//////////////////////////////////////////////////////////////////////////////
// GENERATE NEAREST_NEIGHBORS LIST
size_t **generate_nearest_neighbors_list(size_t dimension, size_t lattice_side) {
	size_t i, j;
	size_t n_sites = pow(lattice_side, dimension);
	size_t **nearest_neighbor = alloc_2d_ull(n_sites, 2 * dimension);

	for (i = 0; i < n_sites; i++) {
		size_t *x = get_coordinates_from_site_index(dimension, lattice_side, i);
		size_t *y = (size_t *)malloc(dimension * sizeof(size_t));
		for (j = 0; j < dimension; j++) y[j] = x[j];
		for (j = 0; j < dimension; j++) {
			if(y[j] == 0) y[j] = lattice_side - 1;
			else y[j] -= 1;
			nearest_neighbor[i][2 * j] = get_site_index_from_coordinates(dimension, lattice_side, y);
			y[j] = x[j];
			if (y[j] == lattice_side - 1) y[j] = 0;
			else y[j] += 1;
			nearest_neighbor[i][2 * j + 1] = get_site_index_from_coordinates(dimension, lattice_side, y);
			y[j] = x[j];
		}
		free(x);
		free(y);
	}

	return nearest_neighbor;
}



//////////////////////////////////////////////////////////////////////////////
// BUILD SNAKE
void build_snake_recursive(size_t dimension, size_t lattice_side, size_t depth, size_t *coordinate, size_t *path, size_t *idx, int forward) {
  size_t i;
  if (depth == dimension) {
    path[(*idx)++] = get_site_index_from_coordinates(dimension, lattice_side, coordinate);
    return;
  }
  if (forward) {
    for (i = 0; i < lattice_side; i++) {
      coordinate[depth] = i;
      build_snake_recursive(dimension, lattice_side, depth + 1, coordinate, path, idx, (i % 2 == 0) ? 1 : 0);
    }
  } else {
    for (i = lattice_side; i-- > 0;) {
      coordinate[depth] = i;
      build_snake_recursive(dimension, lattice_side, depth + 1, coordinate, path, idx, (i % 2 == 0) ? 0 : 1);
    }
  }
}


//////////////////////////////////////////////////////////////////////////////
// INITIALIZE CONFIGURATION
Configuration *initialize_configuration(Parameters *param, size_t **nearest_neighbor) {
  size_t i, j;
  size_t d = param->dimension;
  size_t L = param->lattice_side;
  size_t n = param->n_linear_chains;
  size_t V = 1;
  for (i = 0; i < d; i++) V *= L;

  // Check max number of chains
  size_t n_max = (V % 2 == 0) ? (V / 2 - 1) : ((V - 1) / 2);
  if (n > n_max) {
    printf("initialize_configuration: n_linear_chains must be <= %zu\n", n_max);
    return NULL;
  }

  // Allocate configuration
  Configuration *config = malloc(sizeof(Configuration));
  if (!config) return NULL;
  config->dimension = d;
  config->lattice_side = L;
  config->n_sites = V;
  config->n_linear_chains = n;
  config->connected_site = alloc_2d_ull(V, 2);
  config->endpoint = malloc(2 * n * sizeof(size_t));
  if (!config->connected_site || !config->endpoint) return NULL;

  // Initialize connections
  for (i = 0; i < V; i++) {
    config->connected_site[i][0] = SIZE_MAX;
    config->connected_site[i][1] = SIZE_MAX;
  }

  // Build Hamiltonian snake path
  size_t *path = malloc(V * sizeof(size_t));
  size_t *coordinate = calloc(d, sizeof(size_t));
  size_t idx = 0;
  build_snake_recursive(d, L, 0, coordinate, path, &idx, 1);
  free(coordinate);

  // Connect path
  for (i = 0; i < V - 1; i++) {
    size_t a = path[i];
    size_t b = path[i + 1];
    size_t *csA = config->connected_site[a];
    size_t *csB = config->connected_site[b];
    csA[(csA[0] == SIZE_MAX) ? 0 : 1] = b;
    csB[(csB[0] == SIZE_MAX) ? 0 : 1] = a;
  }

  // Cut (n-1) bonds
  for (j = 1; j < n; j++) {
    size_t cut = (j * V) / n;
    size_t a = path[cut - 1];
    size_t b = path[cut];
    size_t *csA = config->connected_site[a];
    size_t *csB = config->connected_site[b];
    if (csA[0] == b) csA[0] = SIZE_MAX;
    else csA[1] = SIZE_MAX;
    if (csB[0] == a) csB[0] = SIZE_MAX;
    else csB[1] = SIZE_MAX;
  }

  // Collect endpoints
  size_t e = 0;
  for (i = 0; i < V; i++) if (config->connected_site[i][0] == SIZE_MAX || config->connected_site[i][1] == SIZE_MAX) config->endpoint[e++] = i;

  free(path);

  return config;
}


//////////////////////////////////////////////////////////////////////////////
// INITIALIZE CONFIGURATION FROM INPUT FILE
Configuration *copy_configuration_from_config_file(const char *config_filename, Parameters *param) {
	if (param->n_linear_chains >= pow(param->lattice_side, param->dimension - 1)) {
		printf("initialize_configuration: n_linear_chains must be less than %zu\n", pow(param->lattice_side, param->dimension - 1));
		return NULL;
	}

  FILE *file = fopen(config_filename, "r");
  if (!file) {
    printf("Error: Could not open file %s\n", config_filename);
    return NULL;
  }

  Configuration *config = (Configuration *)malloc(sizeof(Configuration));
  if (!config) {
    perror("Failed to allocate memory for the configuration\n");
    fclose(file);
    return NULL;
  }

	config->dimension = param->dimension;
	config->lattice_side = param->lattice_side;
	config->n_sites = pow(config->lattice_side, config->dimension);
	config->n_linear_chains = param->n_linear_chains;
	config->endpoint = (size_t *)malloc(2 * param->n_linear_chains * sizeof(size_t));
	config->connected_site = alloc_2d_ull(config->n_sites, 2);

  size_t i;
	size_t j = 0;
  for (i = 0; i < config->n_sites; i++) {
    if (fscanf(file, "%zu %zu", &config->connected_site[i][0], &config->connected_site[i][1]) == EOF) {
      printf("initialize_configuration_from_config_filename: fatal error while reading input file.\n");
      exit(1);
    }
		if (config->connected_site[i][0] == SIZE_MAX || config->connected_site[i][1] == SIZE_MAX) {
			config->endpoint[j] = i;
			j++;
		}
  }

  return config;
}


//////////////////////////////////////////////////////////////////////////////
// CONFIGURATION DESTROYER
void free_configuration(Configuration *config) {
  if (config) {
    if (config->endpoint) free(config->endpoint);
    if (config->connected_site) free_2d_ull(config->connected_site);
    free(config);
  }
}


//////////////////////////////////////////////////////////////////////////////
// SAMC MONTE CARLO MOVE
void SAMC_MOVE(Configuration *config, size_t **nearest_neighbor, long *idum) {
  size_t site_A, site_B, site_C;

  // Choose reconnecting endpoint (site A)
  size_t e = rand_int(2 * config->n_linear_chains, idum);
  site_A = config->endpoint[e];
  size_t **conn = config->connected_site;
  size_t *csA = conn[site_A];

  // Choose unbonded nearest-neighbour of A (site B)
  site_B = nearest_neighbor[site_A][rand_int(2 * config->dimension, idum)];
  while (site_B == csA[0] || site_B == csA[1]) site_B = nearest_neighbor[site_A][rand_int(2 * config->dimension, idum)];

  // Choose site to disconnect from B (site C)
  size_t *csB = conn[site_B];
  size_t c0 = csB[0];
  size_t c1 = csB[1];
  if (c0 == SIZE_MAX) site_C = c1;
  else if (c1 == SIZE_MAX) site_C = c0;
  else site_C = (ran2(idum) < 0.5) ? c0 : c1; // B is not endpoint
	size_t *csC = conn[site_C];
	if (csC[0] == SIZE_MAX || csC[1] == SIZE_MAX) return; // reject if C is endpoint

  // Update endpoint and connected sites
  config->endpoint[e] = site_C;
  if (csA[0] == SIZE_MAX) csA[0] = site_B;
  else csA[1] = site_B;
  if (csB[0] == site_C) csB[0] = site_A;
  else csB[1] = site_A;
  if (csC[0] == site_B) csC[0] = SIZE_MAX;
  else csC[1] = SIZE_MAX;
}


//////////////////////////////////////////////////////////////////////////////
// SAMC TRAJECTORY
void samc_trajectory(Configuration *config, size_t **nearest_neighbor, Parameters *param, long *idum, bool output) {
  char *system = (char *)malloc(64 * sizeof(char));
  if (param->dimension == 2) snprintf(system, 64, "%zux%zu", param->lattice_side, param->lattice_side);
  else if (param->dimension == 3) snprintf(system, 64, "%zux%zux%zu", param->lattice_side, param->lattice_side, param->lattice_side);
  else if (param->dimension == 4) snprintf(system, 64, "%zux%zux%zux%zu", param->lattice_side, param->lattice_side, param->lattice_side, param->lattice_side);
  else {
    printf("Error: Unsupported lattice param->dimension %zu\n", param->dimension);
    return;
  }
  
  // Trajectory
  size_t i, j;
	clock_t start, now;
	double cpu_time;
	start = clock();

  for (i = 0; i < param->n_dumps; i++) {
    for (j = 0; j < param->n_steps_per_dump; j++) SAMC_MOVE(config, nearest_neighbor, idum); // Evolve n_steps_per_dump times
    if (output) {
      char *filename = (char *)malloc(256 * sizeof(char));
  		snprintf(filename, 256, "configurations/config_%s_%zulinear_%09zu.dat", system, param->n_linear_chains, i);
      FILE *file = fopen(filename, "w");
      if (!file) {
        printf("Error: Could not open file %s for writing.\n", filename);
        return;
      }
      for (j = 0; j < config->n_sites; j++) fprintf(file, "%zu %zu\n", config->connected_site[j][0], config->connected_site[j][1]);
		  free(filename);
  	  fclose(file);
    }
    // Print on terminal
    if ((i + 1) % (param->n_dumps / 10) == 0) {
			now = clock();
			cpu_time = ((double)(now - start)) / CLOCKS_PER_SEC;
			printf("%zu configurations generated in %g seconds\n", param->n_dumps / 10, cpu_time);
			start = now;
		}
  }
  printf("End of trajectory. %zu configurations have been generated\n", param->n_dumps);

	free(system);  
}


//////////////////////////////////////////////////////////////////////////////
// RANDOM INTEGER
size_t rand_int(size_t max_value, long *idum) { return (size_t)(max_value * ran2(idum)); }


//////////////////////////////////////////////////////////////////////////////
// MEMORY ALLOCATION FUNCTION
// Allocate 2D arrays as pointer table (rows) + one contiguous data block, avoiding O(n_sites) small mallocs and drastically reduces fragmentation/overhead
size_t **alloc_2d_ull(size_t rows, size_t cols) {
	if (rows != 0 && cols > SIZE_MAX / rows) return NULL;
	size_t n = rows * cols;
	if (n > SIZE_MAX / sizeof(size_t)) return NULL;
	size_t **a = malloc(rows * sizeof(*a));
	if (!a) return NULL;
	a[0] = malloc(n * sizeof(**a));
	if (!a[0]) { free(a); return NULL; }
	size_t i;
	for (i = 1; i < rows; i++)
	a[i] = a[0] + i * cols;

	return a;
}


// MEMORY CLEANER
void free_2d_ull(size_t **a) {
	if (!a) return;
	free(a[0]);
	free(a);
}


//////////////////////////////////////////////////////////////////////////////
// EXTRACT CHAIN COORDINATES
void print_unwrapped_configuration(Configuration *config) {
	size_t e, i, d;
	int k;
	size_t dim = config->dimension;
	size_t lattice_side = config->lattice_side;
	size_t n_sites = config->n_sites;
	bool *visited = calloc(n_sites, sizeof(bool));
	size_t n_endpoints = 2 * config->n_linear_chains;

	// LINEAR CHAINS
	for (e = 0; e < n_endpoints; e++) {
		size_t start = config->endpoint[e];
		if (visited[start]) continue;
		size_t prev = SIZE_MAX;
		size_t curr = start;
		long *unwrap = calloc(dim, sizeof(long));
		size_t *coord = get_coordinates_from_site_index(dim, lattice_side, curr);
		for (d = 0; d < dim; d++) unwrap[d] = coord[d];
		size_t length = 0;
		printf("linear ");

		// First pass: count length
		{
			size_t tmp_prev = SIZE_MAX;
			size_t tmp_curr = start;
			while (tmp_curr != SIZE_MAX && !visited[tmp_curr]) {
				visited[tmp_curr] = true;
				length++;
				size_t next = SIZE_MAX;
				int k;
				for (k = 0; k < 2; k++) {
      		size_t c = config->connected_site[tmp_curr][k];
					if (c != SIZE_MAX && c != tmp_prev) {
						next = c;
						break;
					}
				}
				tmp_prev = tmp_curr;
				tmp_curr = next;
			}
		}
		printf("%zu\n", length - 1);

		// Reset visited for proper traversal
		curr = start;
		prev = SIZE_MAX;
		for (i = 0; i < length; i++)	{
			printf("%ld", unwrap[0]);
			for (d = 1; d < dim; d++) printf(" %ld", unwrap[d]);
			printf("\n");
			size_t next = SIZE_MAX;
			for (k = 0; k < 2; k++) {
				size_t c = config->connected_site[curr][k];
				if (c != SIZE_MAX && c != prev) {
					next = c;
					break;
				}
			}
			if (next == SIZE_MAX) break;
			size_t *next_coord = get_coordinates_from_site_index(dim, lattice_side, next);
			for (d = 0; d < dim; d++) {
				int step = periodic_step(coord[d], next_coord[d], lattice_side);
				unwrap[d] += step;
			}
			free(coord);
			coord = next_coord;
			prev = curr;
			curr = next;
		}
		free(coord);
		free(unwrap);
	}

	// RINGS
	for (i = 0; i < n_sites; i++) {
		if (visited[i]) continue;
		size_t start = i;
		size_t prev = SIZE_MAX;
		size_t curr = start;
		long *unwrap = calloc(dim, sizeof(long));
		size_t *coord = get_coordinates_from_site_index(dim, lattice_side, curr);
		for (d = 0; d < dim; d++) unwrap[d] = coord[d];
		long total_disp[dim];
		for (d = 0; d < dim; d++) total_disp[d] = 0;
		size_t length = 0;

		// First pass to count
		do {
			visited[curr] = true;
			length++;
			size_t next = (config->connected_site[curr][0] != prev) ? config->connected_site[curr][0] : config->connected_site[curr][1];
			prev = curr;
			curr = next;
		} while (curr != start);
		printf("ring ");

		// Second pass for unwrapping
		prev = SIZE_MAX;
		curr = start;
		bool winding = false;
		printf("%zu\n", length);

		do {
			printf("%ld", unwrap[0]);
			for (d = 1; d < dim; d++) printf(" %ld", unwrap[d]);
			printf("\n");
			size_t next = (config->connected_site[curr][0] != prev) ? config->connected_site[curr][0] : config->connected_site[curr][1];
			size_t *next_coord = get_coordinates_from_site_index(dim, lattice_side, next);
			for (d = 0; d < dim; d++) {
				int step = periodic_step(coord[d], next_coord[d], lattice_side);
				unwrap[d] += step;
				total_disp[d] += step;
			}
			free(coord);
			coord = next_coord;
			prev = curr;
			curr = next;
		} while (curr != start);

		for (d = 0; d < dim; d++) if (total_disp[d] != 0) winding = true;
		if (winding) printf("#winding\n"); // written below coordinates!

		free(coord);
		free(unwrap);
	}

	free(visited);
}


// PERIODIC STEP
static int periodic_step(size_t a, size_t b, size_t lattice_side) {
	if ((a + 1) % lattice_side == b) return +1;
	if ((b + 1) % lattice_side == a) return -1;
	return 0;
}
