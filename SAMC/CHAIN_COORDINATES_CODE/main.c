////////////////////////////////////////////////
// efornasa, October 2025                     //
// Processing configurations: extract chains' //
// coordinates and metric properties.         //
////////////////////////////////////////////////

// main.c

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <math.h>
#include "functions.h"


int main(int argc, char *argv[]) {
  if (argc != 5) {
    fprintf(stderr, "Usage: %s <system> <n_linear_chains> <start_config> <n_configs>\n", argv[0]);
    return EXIT_FAILURE;
  }

  const char *system = argv[1];
	int n_linear_chains = atoi(argv[2]);
  int start_config = atoi(argv[3]);
  int n_configs = atoi(argv[4]);
  
  // Create the lattice
  char *lattice_filename = (char *)malloc(64 * sizeof(char));
	sprintf(lattice_filename, "../lattice_files/lattice_%s.dat", system);
  HC_Lattice *lattice = construct_lattice_from_input_file(lattice_filename);
  if (!lattice) {
    fprintf(stderr, "Failed to create lattice\n");
    return EXIT_FAILURE;
  }
  printf("Lattice built properly: n_sites = %d, n_edges = %d\n\n", lattice->n_sites, lattice->n_edges);

	free(lattice_filename);

/*
  // Open output files
  char *linear_chains_filename = (char *)malloc(128 * sizeof(char));
  sprintf(linear_chains_filename, "metric_linear_%s_%dlinear.dat", system, n_linear_chains);
  FILE *linear_chains_file = fopen(linear_chains_filename, "w");
  if (!linear_chains_file) {
    perror("Error opening output files");
    return EXIT_FAILURE;
  }
	fprintf(linear_chains_file, "#");
  for (i = 0; i < n_linear_chains; i++) fprintf(linear_chains_file, "l_%d Ree_%d Rg_%d ", i, i, i);
	fprintf(linear_chains_file, "n_rings\n");
*/

  // PROCESS CONFIGURATIONS  
  int i, j, k, l;

  for (i = 0; i < n_configs; i++) {
		// Read configuration file
  	char *config_filename = (char *)malloc(128 * sizeof(char));
		sprintf(config_filename, "configurations/config_%s_%06d.dat", system, start_config + i);

		FILE *config_file = fopen(config_filename, "r");
    if (!config_file) {
      perror("Error opening configuration file");
      return EXIT_FAILURE;
    }

		int *end = (int *)malloc(2 * n_linear_chains * sizeof(int));
    char *edge_string = (char *)malloc((lattice->n_edges + 1) * sizeof(char));
    if (fscanf(config_file, "%s", edge_string) == EOF) {
      printf("Fatal error while reading input file.\n");
      exit(1);
    }
		for (j = 0; j < 2 * n_linear_chains; j++) {
      if (fscanf(config_file, "%d", &end[j]) == EOF) {
        printf("Fatal error while reading input file.\n");
        exit(1);
      }
    }
		fclose(config_file);
		free(config_filename);

    // Check if the string length matches the expected lattice->n_edges
    int edge_string_length = strlen(edge_string);
    if (edge_string_length != lattice->n_edges) {
      fprintf(stderr, "Failed to read the edge string. String length is %d instead of %d\n", edge_string_length, lattice->n_edges);
      break;
    }

    // Open coordinates output file and write header
    char *chain_coordinates_filename = (char *)malloc(128 * sizeof(char));
    sprintf(chain_coordinates_filename, "chain_coordinates/chain_coordinates_%s_%06d.dat", system, start_config + i);

    FILE *chain_coordinates_file = fopen(chain_coordinates_filename, "w");
    if (!chain_coordinates_file) {
      perror("Error opening output files");
      return EXIT_FAILURE;
    }
    fprintf(chain_coordinates_file, "#line format: site_index coord[0] ... coord[d-1]\n");

    // Find chains and rings
    Configuration *config = construct_configuration_chains(lattice, edge_string, end, n_linear_chains);
    if (!config) {
      fprintf(stderr, "Failed to construct chains and rings.\n");
      continue;
    }
		free(edge_string);
		free(end);
 
/*   
    // Write out OPEN CHAINS metric properties
		for (j = 0; j < n_linear_chains; j++) fprintf(linear_chains_file, "%d %g %g ", config->linear_chain[j]->length, config->linear_chain[j]->Ree, config->linear_chain[j]->Rg);
		fprintf(linear_chains_file, "%d\n", n_rings);
*/

    // Write out COORDINATES
		for (j = 0; j < n_linear_chains; j++) {
    	fprintf(chain_coordinates_file, "linear %d\n", config->linear_chain[j]->length);
      for (k = 0; k <= config->linear_chain[j]->length; k++) {
        fprintf(chain_coordinates_file, "%d", config->linear_chain[j]->chain_array[k]);
        for (l = 0; l < lattice->dimension; l++) fprintf(chain_coordinates_file, " %d", config->linear_chain[j]->chain_coordinate[k][l]);
        fprintf(chain_coordinates_file, "\n");
      }
		}

    for (j = 0; j < config->n_rings; j++) {
      fprintf(chain_coordinates_file, "ring %d\n", config->ring[j]->length);
      for (k = 0; k < config->ring[j]->length; k++) {
        fprintf(chain_coordinates_file, "%d", config->ring[j]->chain_array[k]);
        for (l = 0; l < lattice->dimension; l++) fprintf(chain_coordinates_file, " %d", config->ring[j]->chain_coordinate[k][l]);
        fprintf(chain_coordinates_file, "\n");
      }
    }
    free_configuration(config);

    fclose(chain_coordinates_file);
    free(chain_coordinates_filename);
  }
  
  // Free memory
  free_lattice(lattice);

  return EXIT_SUCCESS;
}
