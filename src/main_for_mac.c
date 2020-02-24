#include "neighborhood_search_for_mac.h"

int main()
{

	GLfloat(*data)[8] = malloc(sizeof(data[0]) * NPTS);
	CHECK_MALLOC(data);
	GLfloat(*coord)[2] = malloc(sizeof(coord[0]) * NPTS);
	CHECK_MALLOC(coord);

	// Seed the random
	time_t seed = time(NULL);
	printf(" %u \n", seed);
	srand(seed);
	fillData_for_mac(data, coord, NPTS);

	neighborhood_for_mac** nhList = create_neighborhood_for_mac(data, coord, MAX_ITER, 1, NPTS, USE_CELLS, USE_IMPROVED_METHOD, RADIUS_ALGORITHM, SEARCH_NEIGHBORHOOD, USE_VERLET, seed);

	for (int i = 0; i < MAX_ITER; i++)
		neighborhood_delete_for_mac(nhList[i], NPTS);
	free(nhList);
	free(coord);
	free(data);

	show_equivalence_of_major_algorithm_combination_for_mac(0);

	return EXIT_SUCCESS;
}