#include "neighborhood_search.h"

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
	fillData(data, coord, NPTS);

	neighborhood** nhList = create_neighborhood(data, coord, MAX_ITER, 1, NPTS, USE_CELLS, USE_IMPROVED_METHOD, RADIUS_ALGORITHM, SEARCH_NEIGHBORHOOD, USE_VERLET, USE_THREADS, NTHREADS, seed);

	for (int i = 0; i < MAX_ITER; i++)
		neighborhood_delete(nhList[i], NPTS);
	free(nhList);
	free(coord);
	free(data);

	show_equivalence_of_major_algorithm_combination(0);

	return EXIT_SUCCESS;
}