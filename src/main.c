#include "BOV.h"
#include <time.h>
#include <math.h>
#include <pthread.h>
#include <pthread.c>
#include <implement.h>

// see stringification process
#define xstr(s) str(s)
#define str(s) #s
#define M_PI 3.14159265358979323846

// parameters
#define NPTS 25000
#define USE_CELLS 1
#define USE_IMPROVED_METHOD 1
#define RADIUS_ALGORITHM 1
#define USE_VERLET 1
#define USE_THREADS 1
#define nThreads 5
#define MAX_ITER 10
#define KH 0

// malloc verification
#define CHECK_MALLOC(ptr) if((ptr)==NULL) { \
		BOV_ERROR_LOG(BOV_OUT_OF_MEM_ERROR, "Memory allocation failed"); \
		exit(EXIT_FAILURE); }

typedef struct neighbours {
	int index;
	float distance;
	struct neighbours* next;
}neighbours;

typedef struct neighborhood {
	pthread_mutex_t mutex;
	int index;
	int nNeighbours;
	int nPotentialNeighbours;
	neighbours* list;
	neighbours* potential_list;
}neighborhood;

typedef struct node {
	int index;
	struct node* next;
}node;

typedef struct cell {
	int nResident;
	node* ResidentList;
}cell;

typedef struct protected_int {
	int var;
	pthread_mutex_t mutex;
}protected_int;

typedef struct loop_thread_arg {
	cell* cells;
	protected_int* cellCounter;
	GLsizei nPoints;
	int iterations;
	double size;
	neighborhood* nh;
	double kh;
	double L;
	GLfloat(*coord)[2];
	GLfloat(*data)[8];
	int this_thread_number;
	int NTH;
	int use_verlet;
	int use_cells;
}lta;


struct neighbours* neighbours_new(int index, neighborhood* neigh, int i, float d, int potential)
{
	neighbours* new = malloc(sizeof(neighbours));
	CHECK_MALLOC(new);
	new->index = index;
	new->distance = d;
	if (potential) {
		pthread_mutex_lock(&neigh[i].mutex);
		new->next = neigh[i].potential_list;
		neigh[i].potential_list = new;
		neigh[i].nPotentialNeighbours++;
		pthread_mutex_unlock(&neigh[i].mutex);
	}
	else {
		pthread_mutex_lock(&neigh[i].mutex);
		new->next = neigh[i].list;
		neigh[i].list = new;
		neigh[i].nNeighbours++;
		pthread_mutex_unlock(&neigh[i].mutex);
	}
	return new;
}

void neighbours_delete(neighbours* n) {
	if (n) {
		neighbours* temp = n->next;
		free(n);
		neighbours_delete(temp);
	}
}

neighborhood* neighborhood_new(GLsizei n)
{
	neighborhood* neigh = calloc(n, sizeof(neighborhood));
	CHECK_MALLOC(neigh);
	for (int i = 0; i < n; i++) {
		neigh[i].index = i;
		neigh[i].nNeighbours = 0;
		neigh[i].nPotentialNeighbours = 0;
		neigh[i].list = NULL;
		neigh[i].potential_list = NULL;
		neigh[i].mutex = PTHREAD_MUTEX_INITIALIZER;
	}
	return neigh;
}

void neighborhood_delete(neighborhood* nh, GLsizei n) {
	for (int i = 0; i < n; i++) {
		neighbours_delete(nh[i].list);
		neighbours_delete(nh[i].potential_list);
		pthread_mutex_destroy(&nh[i].mutex);
	}
	free(nh);
}

struct node* node_new(cell* c, int i, int index)
{
	node* new = malloc(sizeof(node));
	CHECK_MALLOC(new);
	new->index = index;
	new->next = c[i].ResidentList;
	c[i].ResidentList = new;
	c[i].nResident++;
	return new;
}

void node_delete(node* n) {
	if (n) {
		neighbours* temp = n->next;
		free(n);
		node_delete(temp);
	}
}

cell* cell_new(GLsizei n)
{
	cell* c = calloc(n, sizeof(cell));
	CHECK_MALLOC(c);
	for (int i = 0; i < n; i++) {
		c[i].nResident = 0;
		c[i].ResidentList = NULL;
	}
	return c;
}

void cell_delete(cell* c, GLsizei n) {
	for (int i = 0; i < n; i++)
		node_delete(c[i].ResidentList);
	free(c);
}

void printNeighborhood(neighborhood* nh, GLfloat data[][8], int size) {
	for (int i = 0; i < size; i++) {
		printf("Resident %i : coordinate: %f %f   number of neighbours %i\n", i + 1, data[i][0], data[i][1], nh[i].nNeighbours);
		int j = 1;
		for (neighbours* current = nh[i].list; current; current = current->next)
			printf("   Neighbours %i : %f %f\n", j++, data[current->index][0], data[current->index][1]);
	}
}

void printCell(GLfloat data[][8], cell* c, int size) {
	for (int i = 0; i < size; i++) {
		printf("Cell %i : %i\n", i + 1, c[i].nResident);
		int j = 1;
		for (node* current = c[i].ResidentList; current; current = current->next)
			printf("   Neighbours %i : %f %f\n", j++, data[current->index][0], data[current->index][1]);
	}
}

void protected_inc(protected_int* guard) {
	pthread_mutex_lock(&guard->mutex);
	guard->var++;
	pthread_mutex_unlock(&guard->mutex);
}

int protected_get(protected_int* guard) {
	pthread_mutex_lock(&guard->mutex);
	int value = guard->var;
	pthread_mutex_unlock(&guard->mutex);
	return value;
}

int protected_get_inc(protected_int* guard) {
	pthread_mutex_lock(&guard->mutex);
	int value = guard->var;
	guard->var++;
	pthread_mutex_unlock(&guard->mutex);
	return value;
}

void protected_reset(protected_int* guard) {
	pthread_mutex_lock(&guard->mutex);
	guard->var = 0;
	pthread_mutex_unlock(&guard->mutex);
}

double compute_kh(int nPoints, int RA) {
	double target = 21.0 / nPoints;
	if (nPoints < 21)
		return sqrt(2.0);
	else if (!RA)
		return sqrt(2.0) * target;
	else if (target <= M_PI / 4)
		return sqrt(4 * target / M_PI);
	double tolerance = 0.0000000001;
	double kh_min = 1.0;
	double kh_max = sqrt(2.0);
	while (fabs(kh_max - kh_min) >= tolerance) {
		kh_max = kh_min;
		kh_min = sqrt((target - sin(acos(1.0 / kh_min)) * kh_min) / ((M_PI / 4 - acos(1.0 / kh_min))));
	}
	return kh_min;
}

void* loop_thread(void* args) {
	lta* loop = (lta*)args;
	cell* cellArray = loop->cells;
	protected_int* cellCounter = loop->cellCounter;
	int size = ceil(loop->size);
	int this_thread_number = loop->this_thread_number;
	neighborhood* nh = loop->nh;
	double kh = loop->kh;
	double L = loop->L;
	GLsizei nPoints = loop->nPoints;
	int iterations = loop->iterations;
	int use_verlet = loop->use_verlet;
	int use_cells = loop->use_cells;
	int NTH = loop->NTH;
	unsigned long frameCount = 0;
	int i;
	if (use_cells)
		i = 0;
	else
		i = (int)(this_thread_number * ((double)nPoints) / NTH);
	int j = 1;
	int nValid = 0;
	int nInvalid = 0;
	int nVerlet = 0;
	int checking_cell_number = -1;
	int this_cell_number = -1;
	int checked_cells = 0;
	cell checking_cell;
	node checking_node;
	neighbours checking_neighbours;
	cell this_cell;
	node this_node;
	int i_check = i - 1;
	int j_check = j - 1;
	int f_check = 1;
	int potential_or_not = 0;
	while (((!use_cells && i < ((int)((this_thread_number + 1) * ((double)nPoints) / NTH))) || (use_cells && i < nPoints)) && this_cell_number < size * size) {
		if (i != i_check) {
			nValid = 0;
			nInvalid = 0;
			nVerlet = 0;
			if (use_cells) {
				if (i == 0) {
					this_cell_number = protected_get_inc(cellCounter);
					while (!cellArray[this_cell_number].nResident && this_cell_number < size * size)
						this_cell_number = protected_get_inc(cellCounter);
					if (this_cell_number < size * size) {
						this_cell = cellArray[this_cell_number];
						this_node = this_cell.ResidentList[0];
					}
				}
				else {
					if (this_node.next)
						this_node = *this_node.next;
					else {
						this_cell_number = protected_get_inc(cellCounter);
						while (this_cell_number < size * size && !cellArray[this_cell_number].nResident && !cellArray[this_cell_number].ResidentList)
							this_cell_number = protected_get_inc(cellCounter);
						if (this_cell_number < size * size) {
							this_cell = cellArray[this_cell_number];
							this_node = this_cell.ResidentList[0];
						}
					}
				}
				if (use_verlet && iterations) {
					if (nh[this_node.index].list) {
						potential_or_not = 0;
						checking_neighbours = nh[this_node.index].list[0];
					}
					else if (nh[this_node.index].potential_list) {
						potential_or_not = 1;
						checking_neighbours = nh[this_node.index].potential_list[0];
					}
					else
						checking_cell_number = -1;
				}
				else {
					checked_cells = 0;
					if (USE_IMPROVED_METHOD && this_cell_number < size * size) {
						if (this_node.next) {
							checking_cell_number = this_cell_number;
							checking_cell = cellArray[checking_cell_number];
							checking_node = *this_node.next;
						}
						else {
							checking_cell_number = find_next_cell(this_cell_number, ++checked_cells, size, USE_IMPROVED_METHOD);
							while (checking_cell_number != -1 && !cellArray[checking_cell_number].nResident && !cellArray[checking_cell_number].ResidentList)
								checking_cell_number = find_next_cell(this_cell_number, ++checked_cells, size, USE_IMPROVED_METHOD);
							if (checking_cell_number != -1) {
								checking_cell = cellArray[checking_cell_number];
								checking_node = checking_cell.ResidentList[0];
							}
						}
					}
					else if (this_cell_number < size * size) {
						checking_cell_number = find_next_cell(this_cell_number, checked_cells, size, USE_IMPROVED_METHOD);
						while (checking_cell_number != -1 && !cellArray[checking_cell_number].nResident && !cellArray[checking_cell_number].ResidentList)
							checking_cell_number = find_next_cell(this_cell_number, ++checked_cells, size, USE_IMPROVED_METHOD);
						if (checking_cell_number != -1) {
							checking_cell = cellArray[checking_cell_number];
							checking_node = checking_cell.ResidentList[0];
						}
					}
				}
			}
			else if (use_verlet && iterations) {
				if (nh[i].list) {
					potential_or_not = 0;
					checking_neighbours = nh[i].list[0];
				}
				else if (nh[i].potential_list) {
					potential_or_not = 1;
					checking_neighbours = nh[i].potential_list[0];
				}
				else
					checking_cell_number = -1;
			}
			if (USE_IMPROVED_METHOD) {
				j_check = i;
				j = i + 1;
			}
		}
		i_check = i;
		if ((use_cells || (use_verlet && iterations)) && checking_cell_number == -1) {
			i++;
			j_check = j;
		}
		if (j != j_check) {
			int index_i, index_j;
			if (use_verlet && iterations) {
				index_j = checking_neighbours.index;
				index_i = this_node.index;
			}
			else if (use_cells) {
				index_j = checking_node.index;
				index_i = this_node.index;
			}
			else {
				index_j = j;
				index_i = i;
			}
			float distance = (pow((double)loop->coord[index_j][0] - (double)loop->coord[index_i][0], 2) + pow((double)loop->coord[index_j][1] - (double)loop->coord[index_i][1], 2));
			if (distance <= pow(kh, 2)) {
				neighbours_new(loop->coord[index_j], nh, index_i, distance, 0);
				if (USE_IMPROVED_METHOD)
					neighbours_new(loop->coord[index_i], nh, index_j, distance, 0);
			}
			else if (use_verlet && distance <= pow(kh + L, 2)) {
				neighbours_new(loop->coord[index_j], nh, index_i, distance, 1);
				if (USE_IMPROVED_METHOD)
					neighbours_new(loop->coord[index_i], nh, index_j, distance, 1);
			}
			if (use_verlet && iterations) {
				if (checking_neighbours.next)
					checking_neighbours = *(checking_neighbours.next);
				else if (!potential_or_not)
					checking_neighbours = nh[this_node.index].potential_list[0];
				else
					i++;
			}
			else if (use_cells)
				if (checking_node.next)
					checking_node = *(checking_node.next);
				else {
					checking_cell_number = find_next_cell(this_cell_number, ++checked_cells, size, USE_IMPROVED_METHOD);
					while (checking_cell_number != -1 && !cellArray[checking_cell_number].nResident && !cellArray[checking_cell_number].ResidentList)
						checking_cell_number = find_next_cell(this_cell_number, ++checked_cells, size, USE_IMPROVED_METHOD);
					if (checking_cell_number != -1) {
						checking_cell = cellArray[checking_cell_number];
						checking_node = cellArray[checking_cell_number].ResidentList[0];
					}
					else
						i++;
				}
			else
				if (j == nPoints - 1 || (i == nPoints - 1 && (j == nPoints - 2 || USE_IMPROVED_METHOD))) {
					i++;
				}
		}
		j_check = j;
		if (USE_IMPROVED_METHOD) {
			j++;
		}
		else if (!USE_IMPROVED_METHOD)
			j = (frameCount) % (nPoints - 1) + (int)(i <= ((frameCount) % (nPoints - 1)));
		frameCount++;
	}
}

int find_next_cell(int this_cell, int checked_cells, int size, int improve) {
	if (this_cell == 0) {
		if (checked_cells != 4)
			return this_cell + checked_cells / 2 * size + checked_cells % 2;
	}
	else if (this_cell == size - 1) {
		checked_cells += improve;
		if (checked_cells != 4)
			return this_cell - 1 + checked_cells / 2 * size + checked_cells % 2;
	}
	else if (this_cell == (size - 1) * size) {
		checked_cells += 2 * improve;
		if (checked_cells != 4)
			return this_cell - size + checked_cells / 2 * size + checked_cells % 2;
	}
	else if (this_cell == size * size - 1) {
		checked_cells += 3 * improve;
		if (checked_cells < 4)
			return this_cell - size - 1 + checked_cells / 2 * size + checked_cells % 2;
	}
	else if (this_cell / size == 0) {
		checked_cells += improve;
		if (checked_cells != 6)
			return this_cell - 1 + checked_cells / 3 * size + checked_cells % 3;
	}
	else if (this_cell / size == size - 1) {
		checked_cells += 4 * improve;
		if (checked_cells != 6)
			return this_cell - size - 1 + checked_cells / 3 * size + checked_cells % 3;
	}
	else if (this_cell % size == 0) {
		checked_cells += 2 * improve;
		if (checked_cells != 6)
			return this_cell - size + checked_cells / 2 * size + checked_cells % 2;
	}
	else if (this_cell % size == size - 1) {
		checked_cells += 3 * improve;
		if (checked_cells != 6)
			return this_cell - size - 1 + checked_cells / 2 * size + checked_cells % 2;
	}
	else {
		checked_cells += 4 * improve;
		if (checked_cells != 9)
			return this_cell - size - 1 + checked_cells / 3 * size + checked_cells % 3;
	}
	return -1;
}

/*Fills the data vector*/
static void fillData(GLfloat(*data)[8], GLfloat(*coord)[2])
{
	float rmax = 100.0 * sqrtf(2.0f);
	for (int i = 0; i < NPTS; i++) {
		data[i][0] = rand() * 200.0 / RAND_MAX - 100.0; // x (rand between -100 and 100)
		data[i][1] = rand() * 200.0 / RAND_MAX - 100.0; // y (rand between -100 and 100)
		coord[i][0] = data[i][0];
		coord[i][1] = data[i][1];
		float r = sqrt(data[i][0] * data[i][0] + data[i][1] * data[i][1]);
		//data[i][2] = 0; //speed x (not used by default visualization)
		//data[i][3] = 0; // speed y (not used by default visualization)
		data[i][2] = rand() * 2.0 / RAND_MAX - 1.0; //Random starting speed
		data[i][3] = rand() * 2.0 / RAND_MAX - 1.0; //Random starting speed
		data[i][4] = 0.0;
		data[i][5] = 0.0;
		data[i][6] = 0.0;
		data[i][7] = 0.8f; // transparency
	}
}

/*Changes the particle velocities randomly and updates the positions based, we assume elastic collisions with boundaries:
	-timestep: time intervals at which these are updated
	-xmin,xmax,ymin,ymax: boundaries of the domain
	-maxspeed: the maximum speed that can be reached by the particles
*/
void bouncyrandomupdate(GLfloat data[][8], GLfloat coord[][2], double timestep, double xmin, double xmax, double ymin, double ymax, double maxspeed) {
	for (int i = 0; i < NPTS; i++) {
		double speed = sqrtf(data[i][2] * data[i][2] + data[i][3] * data[i][3]);
		data[i][0] += data[i][2] * timestep;
		data[i][1] += data[i][3] * timestep;
		coord[i][0] = data[i][0];
		coord[i][1] = data[i][1];
		data[i][2] += ((double)rand() / RAND_MAX - 0.5) * (0.05 * maxspeed) * timestep;
		data[i][3] += ((double)rand() / RAND_MAX - 0.5) * (0.05 * maxspeed) * timestep;

		if (speed > maxspeed) {//Slows down if speed too high
			data[i][2] = data[i][2] * 0.9;
			data[i][3] = data[i][3] * 0.9;
		}

		/*This next part of the code handles the cases where a particle bounces of the wall*/

		//Particle is too high in x
		if (data[i][0] >= xmax) {
			data[i][0] -= 2 * (data[i][0] - xmax);
			data[i][2] = -data[i][2];
		}
		//Particle is too low in x
		if (data[i][0] <= xmin) {
			data[i][0] -= 2 * (data[i][0] - xmin);
			data[i][2] = -data[i][2];
		}
		//Particle is too high in y
		if (data[i][1] >= ymax) {
			data[i][1] -= 2 * (data[i][1] - ymax);
			data[i][3] = -data[i][3];
		}
		//Particle is too low in y
		if (data[i][1] <= ymin) {
			data[i][1] -= 2 * (data[i][1] - ymin);
			data[i][3] = -data[i][3];
		}
	}
}

int compute_optimal_verlet(double timestep, double maxspeed, double kh) {
	//Boils down to solving a cubic function
	double a = -M_PI * 4 * timestep * timestep * maxspeed * maxspeed;

	double b = -4 * (M_PI * timestep * kh * maxspeed - 9 * maxspeed * maxspeed * timestep * timestep);
	double c = kh * kh * (9 - M_PI) - 36 * timestep * kh * maxspeed;
	double d = -9 * kh * kh;
	double adev = a * 3;
	double bdev = b * 2;
	double cdev = c;
	double discri = bdev * bdev - 4 * adev * cdev;
	if (discri <= 0) {
		return -1;
	}
	else {
		double first_ans = (-bdev + sqrt(discri)) / (2 * adev);
		double second_ans = (-bdev - sqrt(discri)) / (2 * adev);
		double first_y = a * first_ans * first_ans * first_ans + b * first_ans * first_ans + c * first_ans + d;
		double second_y = a * second_ans * second_ans * second_ans + b * second_ans * second_ans + c * second_ans + d;
		double max_ans;
		if (first_y > second_y) {
			max_ans = first_ans;
		}
		else {
			max_ans = second_ans;
		}
		if (max_ans < 2) { return -1; }
		else { return max(ceil(max_ans), floor(max_ans)); }
	}
}

int main()
{
	//Defining the domain
	int xmax = 100;
	int xmin = -100;
	int ymax = 100;
	int ymin = -100;

	const GLsizei nPoints = NPTS;
	double points_width = 0.02;

	double timestep = 0.5;
	double maxspeed = 1;
	int optimal_verlet_steps;

	// Seed the random 
	srand(time(NULL));

	pthread_t* threads = NULL;
	int use_thread = USE_THREADS && nThreads >= 1;
	int NTH = nThreads;
	if (use_thread) {
		if (nPoints < nThreads)
			NTH = nPoints;
	}
	else {
		use_thread = 1;
		NTH = 1;
	}
	threads = malloc(NTH * sizeof(pthread_t));
	double kh = KH;
	if(kh == 0)
		kh = compute_kh(nPoints, RADIUS_ALGORITHM) * (xmax - xmin);
	double L = 0.0;
	int use_verlet = USE_VERLET;
	if (USE_VERLET == 1) {
		optimal_verlet_steps = compute_optimal_verlet(timestep, maxspeed, kh);
		if (optimal_verlet_steps == -1) { use_verlet = 0; }
		else { L = optimal_verlet_steps * timestep * maxspeed; }
	}

	GLfloat(*data)[8] = malloc(sizeof(data[0]) * NPTS);
	GLfloat(*coord)[2] = malloc(sizeof(coord[0]) * NPTS);
	fillData(data, coord);

	// create the gridded data structure
	double width = (xmax - xmin) / (kh + L);
	double height = (ymax - ymin) / (kh + L);

	cell* cellArray = NULL;
	double size = (xmax - xmin) / (kh + L);
	double cellLength = (xmax - xmin) / (size);
	int use_cells = USE_CELLS && (kh + L) < (xmax - xmin) / 3.0;

	protected_int* cellCounter = NULL;
	if (use_cells) {
		cellCounter = malloc(sizeof(protected_int));
		cellCounter->mutex = PTHREAD_MUTEX_INITIALIZER;
		cellCounter->var = 0;
	}
	int iterations = 0;
	while (iterations < MAX_ITER) {
		neighborhood* nh = neighborhood_new(nPoints);

		if (use_cells) {
			protected_reset(cellCounter);
			cellArray = cell_new(ceil(size) * ceil(size));
			for (int i = 0; i < nPoints; i++) {
				int cellNumber = ((int)((coord[i][1] - xmin) / (xmax - xmin) * size) * ceil(size) + (int)((coord[i][0] - xmin) / (xmax - xmin) * size));
				if (data[i][1] == xmax)
					cellNumber -= size;
				if (data[i][0] == xmax)
					cellNumber -= 1;
				node_new(cellArray, cellNumber, i);
			}
		}
		for (int i = 0; i < NTH; i++) {
			lta* loop = malloc(sizeof(lta));
			loop->cells = cellArray;
			loop->use_cells = use_cells;
			loop->use_verlet = use_verlet;
			loop->size = size;
			loop->nh = nh;
			loop->kh = kh;
			loop->L = L;
			loop->nPoints = nPoints;
			loop->coord = coord;
			loop->data = data;
			loop->NTH = NTH;
			loop->this_thread_number = i;
			loop->cellCounter = cellCounter;
			if (use_verlet)
				loop->iterations = iterations % optimal_verlet_steps;
			else
				loop->iterations = 0;
			pthread_create(&threads[i], NULL, loop_thread, loop);
		}
		for (int i = 0; i < NTH; i++) {
			pthread_join(threads[i], NULL);
		}
		// Moves the particles around
		if (use_cells)
			cell_delete(cellArray, ceil(size) * ceil(size));
		neighborhood_delete(nh, NPTS);
		bouncyrandomupdate(data, coord, 1, -100, 100, -100, 100, 2);
		iterations++;
	}

	if (use_cells) {
		pthread_mutex_destroy(&(cellCounter->mutex));
		free(cellCounter);
	}

	free(data);

	return EXIT_SUCCESS;
}