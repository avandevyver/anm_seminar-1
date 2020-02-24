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
#define NPTS 25
#define USE_CELLS 1
#define USE_IMPROVED_METHOD 1
#define RADIUS_ALGORITHM 1
#define DRAWING 1
#define USE_VERLET 1
#define SEARCH_NEIGHBORHOOD 1
#define USE_THREADS 0
#define nThreads 5

// malloc verification
#define CHECK_MALLOC(ptr) if((ptr)==NULL) { \
		BOV_ERROR_LOG(BOV_OUT_OF_MEM_ERROR, "Memory allocation failed"); \
		exit(EXIT_FAILURE); }

// for v, a floating point value between 0 and 1, this function fills color with
// the improved jet colormap color corresponding to v
static void colormap(float v, float color[3])
{
	float v1 = 3.5*(v-0.7);
	float v2 = 1.25*v;
	float v3 = fminf(0.5,v)*2.0;

	color[0] = -v1*v1+1.0f;
	color[1] = 6.0f*v2*v2*(1.0f-v2);
	color[2] = 5.5f*v3*(1.0f-v3)*(1.0f-v3);

	// alternative: classical jet colormap
	// color[0] = 1.5 - 4.0 * fabs(v - 0.75);
	// color[1] = 1.5 - 4.0 * fabs(v - 0.5 );
	// color[2] = 1.5 - 4.0 * fabs(v - 0.25);
}

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

struct neighbours* neighbours_new(int index, neighborhood* neigh, int i, float d, int potential)
{
	neighbours* new = malloc(sizeof(neighbours));
	CHECK_MALLOC(new);
	new->index = index;
	new->distance = d;
	if(potential){
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

struct node* node_new( cell* c, int i, int index)
{
	node* new = malloc(sizeof(node));
	CHECK_MALLOC(new);
	new->index = index;
	new->next = c[i].ResidentList;
	c[i].ResidentList = new;
	c[i].nResident++;
	return new;
}

typedef struct random_points_computation_arg {
	GLfloat* coord[2];
	GLsizei n;
	double scale;
}rpca;

typedef struct protected_int {
	int var;
	pthread_mutex_t mutex;
}protected_int;

typedef struct loop_thread_arg {
	bov_window_t* window;
	GLuint** tabs;
	bov_order_t** orders;
	bov_points_t** points;
	cell* cells;
	protected_int* cellCounter;
	protected_int* finished_threads;
	GLsizei nPoints;
	double scale;
	int div;
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
	guard->var=0;
	pthread_mutex_unlock(&guard->mutex);
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

void printNeighborhood(neighborhood* nh,GLfloat data[][8], int size) {
	for (int i = 0; i < size; i++) {
		printf("Resident %i : coordinate: %f %f   number of neighbours %i\n", i + 1, data[i][0], data[i][1], nh[i].nNeighbours);
		int j = 1;
		for (neighbours* current = nh[i].list; current; current = current->next)
			printf("   Neighbours %i : %f %f\n", j++, data[current->index][0], data[current->index][1]);
	}
}

void printCell(GLfloat data[][8],cell* c, int size) {
	for (int i = 0; i < size; i++) {
		printf("Cell %i : %i\n", i + 1, c[i].nResident);
		int j = 1;
		for (node* current = c[i].ResidentList; current; current = current->next)
			printf("   Neighbours %i : %f %f\n", j++, data[current->index][0], data[current->index][1]);
	}
}

void* random_uniform_points(void* args)
{
	rpca* arguments = (rpca*)args;
	double scale = arguments->scale;
	int n = arguments->n;
	GLfloat** coord = arguments->coord;
	GLfloat min[2] = { -1.0 * scale, -1.0 * scale };
	GLfloat max[2] = { 1.0 * scale, 1.0 * scale };
	for (GLsizei i = 0; i < n; i++) {
		coord[i][0] = (max[0] - min[0]) * rand() / RAND_MAX + min[0];
		coord[i][1] = (max[1] - min[1]) * rand() / RAND_MAX + min[1];
	}
}

void random_uniform_point(GLfloat coord[][2], double scale, int n)
{
	GLfloat min[2] = { -1.0 * scale, -1.0 * scale };
	GLfloat max[2] = { 1.0 * scale, 1.0 * scale };
	for (GLsizei i = 0; i < n; i++) {
		coord[i][0] = (max[0] - min[0]) * rand() / RAND_MAX + min[0];
		coord[i][1] = (max[1] - min[1]) * rand() / RAND_MAX + min[1];
	}
}

void nice_colormap(float color[4], float a)
{
	color[0] = sin(0.33 * a) * 0.3 + 0.7;
	color[1] = sin(0.23 * a + 2.0) * 0.3 + 0.7;
	color[2] = sin(0.17 * a + 5.0) * 0.3 + 0.6;
	color[3] = color[3];
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

bov_points_t* find_grid_points(double kh, double xmax) {
	double size = 2 * xmax / kh;
	float(*pts)[2] = malloc((ceil(size)-1) * 4 * sizeof(pts[0]));
	CHECK_MALLOC(pts);
	for (int i = 1; i < size; i++) {
		pts[(i - 1) * 4][0] = xmax * (2.0 * i / (size) - 1.0);
		pts[(i - 1) * 4][1] = xmax * (-1.0);
		pts[(i - 1) * 4 + 1][0] = xmax * (2.0 * i / (size) - 1.0);
		pts[(i - 1) * 4 + 1][1] = xmax * (1.0);
		pts[(i - 1) * 4 + 2][0] = xmax * (-1.0);
		pts[(i - 1) * 4 + 2][1] = xmax * (2.0 * i / (size) - 1.0);
		pts[(i - 1) * 4 + 3][0] = xmax * (1.0);
		pts[(i - 1) * 4 + 3][1] = xmax * (2.0 * i / (size) - 1.0);
	}
	bov_points_t* points = bov_points_new(pts, (ceil(size) - 1) * 4, GL_STATIC_DRAW);
	return points;
}

void* loop_thread(void* args) {
	lta* loop = (lta*)args;
	cell* cellArray = loop->cells;
	protected_int* cellCounter = loop->cellCounter;
	protected_int* finished_threads = loop->finished_threads;
	int div = loop->div;
	int size = ceil(loop->size);
	int this_thread_number = loop->this_thread_number;
	neighborhood* nh = loop->nh;
	double kh = loop->kh;
	double L = loop->L;
	double scale = loop->scale;
	GLsizei nPoints = loop->nPoints;
	int iterations = loop->iterations;
	int use_verlet = loop->use_verlet;
	int use_cells = loop->use_cells;
	int NTH = loop->NTH;
	double points_width = 0.02;
	char name[32];
	snprintf(name, sizeof(name), "Neighborhood: thread %i", this_thread_number);
	bov_window_t* window=NULL;
	if (DRAWING) {
		window = bov_window_new(400, 400, name);
		bov_window_set_color(window, (GLfloat[]) { 0.9f, 0.85f, 0.8f, 1.0f });
	}
	float test_point_color[4] = { 0.0, 0.0, 1.0, 1.0 };
	float circle_color[4] = { 0.0, 0.0, 1.0, .2 };
	float verlet_circle_color[4] = { 1.0, 0.0, 1.0, .2 };
	float verlet_color[4] = { 1.0, 0.0, 1.0, 1.0 };
	float valid_color[4] = { 0.0, 1.0, 0.0, 1.0 };
	float invalid_color[4] = { 1.0, 0.0, 0.0, 1.0 };
	float other_color[4] = { 0.0, 1.0, 1.0, 1.0 };
	float finished_color[4] = { 1.0, 1.0, 1.0, 1.0 };
	GLuint* tabValid = NULL;
	GLuint* tabInvalid = NULL;
	GLuint* tabFinished = NULL;
	GLuint* tabVerlet = NULL;
	bov_order_t* valid_order = NULL;
	bov_order_t* invalid_order = NULL;
	bov_order_t* finished_order = NULL;
	bov_order_t* verlet_order = NULL;
	bov_points_t* unchecked_points = NULL;
	bov_points_t* this_point = NULL;
	bov_points_t* circle_point = NULL;
	bov_points_t* verlet_circle_point = NULL;
	bov_points_t* checking_point = NULL;
	bov_points_t* valid_points = NULL;
	bov_points_t* invalid_points = NULL;
	bov_points_t* finished_points = NULL;
	bov_points_t* verlet_points = NULL;
	bov_points_t* gridPoint = NULL;
	bov_points_t* pointset = NULL;
	if (DRAWING) {
		tabValid = calloc(nPoints, sizeof(GLuint));
		tabInvalid = calloc(nPoints, sizeof(GLuint));
		tabFinished = calloc(nPoints, sizeof(GLuint));
		valid_order = bov_order_new(tabValid, nPoints + 1, GL_DYNAMIC_DRAW);
		invalid_order = bov_order_new(tabInvalid, nPoints + 1, GL_DYNAMIC_DRAW);
		finished_order = bov_order_new(tabFinished, nPoints + 1, GL_DYNAMIC_DRAW);
		unchecked_points = bov_particles_new(loop->data, nPoints, GL_STATIC_DRAW);
		this_point = bov_points_new(loop->coord, 1, GL_DYNAMIC_DRAW);
		if (use_verlet) {
			tabVerlet = calloc(nPoints, sizeof(GLuint));
			verlet_circle_point = bov_points_new(loop->coord, 1, GL_DYNAMIC_DRAW);
			verlet_order = bov_order_new(tabVerlet, nPoints + 1, GL_DYNAMIC_DRAW);
			verlet_points = bov_points_new(loop->coord, nPoints, GL_STATIC_DRAW);
		}
		circle_point = bov_points_new(loop->coord, 1, GL_DYNAMIC_DRAW);
		checking_point = bov_points_new(loop->coord, 0, GL_DYNAMIC_DRAW);
		valid_points = bov_points_new(loop->coord, nPoints, GL_STATIC_DRAW);
		invalid_points = bov_points_new(loop->coord, nPoints, GL_STATIC_DRAW);
		if(USE_IMPROVED_METHOD)
			finished_points = bov_points_new(loop->coord, nPoints, GL_STATIC_DRAW);
		pointset = bov_points_new((float[4][2]) {
			{-100, 100},
			{ 100, 100 },
			{ 100,-100 },
			{ -100,-100 }
		}, 4, GL_STATIC_DRAW);
		bov_points_param_t lineParams = {
			.fillColor = {1.0, 1.0, 1.0, 1.0},
			.scale = {1.0, 1.0},
			.width = 0.005
		};
		if (use_cells) {
			gridPoint = find_grid_points(kh+L, 100);
			bov_points_set_param(gridPoint, lineParams);
			bov_points_scale(gridPoint, (GLfloat[2]) { scale, scale });
			bov_points_set_pos(gridPoint, (GLfloat[2]) { 0.0, -0.1 });
		}
		if (use_verlet) {
			bov_points_set_color(verlet_circle_point, verlet_circle_color);
			bov_points_set_width(verlet_circle_point, (kh + L) * scale);
			bov_points_scale(verlet_circle_point, (GLfloat[2]) { scale, scale });
			bov_points_set_pos(verlet_circle_point, (GLfloat[2]) { 0.0, -0.1 });
			bov_points_set_color(verlet_points, verlet_color);
			bov_points_set_width(verlet_points, points_width);
			bov_points_set_outline_width(verlet_points, 0.0025);
			bov_points_scale(verlet_points, (GLfloat[2]) { scale, scale });
			bov_points_set_pos(verlet_points, (GLfloat[2]) { 0.0, -0.1 });
		}
		bov_points_set_param(pointset, lineParams);
		bov_points_scale(pointset, (GLfloat[2]) { scale, scale });
		bov_points_set_pos(pointset, (GLfloat[2]) { 0.0, -0.1 });
		bov_points_set_outline_color(unchecked_points, (GLfloat[4]) { 0.3, 0.12, 0.0, 0.25 });
		bov_points_set_color(this_point, test_point_color);
		bov_points_set_color(checking_point, other_color);
		bov_points_set_color(circle_point, circle_color);
		bov_points_set_color(valid_points, valid_color);
		bov_points_set_color(invalid_points, invalid_color);
		if (USE_IMPROVED_METHOD)
		bov_points_set_color(finished_points, finished_color);
		bov_points_set_outline_width(circle_point, 0.);
		bov_points_set_width(circle_point, kh* scale);
		bov_points_set_width(this_point, points_width);
		bov_points_set_width(checking_point, points_width);
		bov_points_set_width(unchecked_points, points_width);
		bov_points_set_width(valid_points, points_width);
		bov_points_set_width(invalid_points, points_width);
		if (USE_IMPROVED_METHOD)
		bov_points_set_width(finished_points, points_width);
		bov_points_set_outline_width(this_point, 0.0025);
		bov_points_set_outline_width(checking_point, 0.0025);
		bov_points_set_outline_width(unchecked_points, 0.0025);
		bov_points_set_outline_width(valid_points, 0.0025);
		bov_points_set_outline_width(invalid_points, 0.0025);
		if (USE_IMPROVED_METHOD)
		bov_points_set_outline_width(finished_points, 0.0025);
		bov_points_scale(this_point, (GLfloat[2]) { scale, scale });
		bov_points_set_pos(this_point, (GLfloat[2]) { 0.0, -0.1 });
		bov_points_scale(checking_point, (GLfloat[2]) { scale, scale });
		bov_points_set_pos(checking_point, (GLfloat[2]) { 0.0, -0.1 });
		bov_points_scale(unchecked_points, (GLfloat[2]) { scale, scale });
		bov_points_set_pos(unchecked_points, (GLfloat[2]) { 0.0, -0.1 });
		bov_points_scale(valid_points, (GLfloat[2]) { scale, scale });
		bov_points_set_pos(valid_points, (GLfloat[2]) { 0.0, -0.1 });
		bov_points_scale(invalid_points, (GLfloat[2]) { scale, scale });
		bov_points_set_pos(invalid_points, (GLfloat[2]) { 0.0, -0.1 });
		if (USE_IMPROVED_METHOD) {
			bov_points_scale(finished_points, (GLfloat[2]) { scale, scale });
			bov_points_set_pos(finished_points, (GLfloat[2]) { 0.0, -0.1 });
		}
		bov_points_scale(circle_point, (GLfloat[2]) { scale, scale });
		bov_points_set_pos(circle_point, (GLfloat[2]) { 0.0, -0.1 });
	}
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
	if (!DRAWING)
		div = 1;
	while ((!DRAWING || !bov_window_should_close(window)) && i < ((int)((this_thread_number + 1) * ((double)nPoints) / NTH)) && this_cell_number < size * size) {
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
				if(use_verlet && iterations){
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
				if (DRAWING && this_cell_number < size * size) {
					bov_points_update(this_point, loop->coord[this_node.index], 1);
					bov_points_update(circle_point, loop->coord[this_node.index], 1);
					if(use_verlet)
						bov_points_update(verlet_circle_point, loop->coord[this_node.index], 1);
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
				if (DRAWING) {
					bov_points_update(this_point, loop->coord[i], 1);
					bov_points_update(circle_point, loop->coord[i], 1);
					if (use_verlet)
						bov_points_update(verlet_circle_point, loop->coord[i], 1);
				}
			}
			else if (DRAWING) {
				bov_points_update(this_point, loop->coord[i], 1);
				bov_points_update(circle_point, loop->coord[i], 1);
				if (use_verlet)
					bov_points_update(verlet_circle_point, loop->coord[i], 1);
			}
			if (USE_IMPROVED_METHOD) {
				if (DRAWING) {
					if (use_cells)
						tabFinished[i] = this_node.index;
					else
						tabFinished[i-(int)(this_thread_number * ((double)nPoints) / NTH)] = i;
					finished_order = bov_order_partial_update(finished_order, tabFinished, 0, i, 0);
				}
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
			if (DRAWING)
				bov_points_update(checking_point, loop->coord[index_j], 1);
			float distance = (pow((double)loop->coord[index_j][0] - (double)loop->coord[index_i][0], 2) + pow((double)loop->coord[index_j][1] - (double)loop->coord[index_i][1], 2));
			if (distance <= pow(kh, 2)) {
				if (DRAWING) {
					tabValid[nValid] = index_j;
					nValid++;
					bov_order_partial_update(valid_order, tabValid, 0, nValid, 0);
				}
				neighbours_new(loop->coord[index_j], nh, index_i, distance,0);
				if (USE_IMPROVED_METHOD)
					neighbours_new(loop->coord[index_i], nh, index_j, distance,0);
			}
			else if (use_verlet && distance <= pow(kh+L, 2)) {
				if (DRAWING) {
					tabVerlet[nVerlet] = index_j;
					nVerlet++;
					bov_order_partial_update(verlet_order, tabVerlet, 0, nVerlet, 0);
				}
				neighbours_new(loop->coord[index_j], nh, index_i, distance,1);
				if (USE_IMPROVED_METHOD)
					neighbours_new(loop->coord[index_i], nh, index_j, distance,1);
			}
			else if (DRAWING) {
				tabInvalid[nInvalid] = index_j;
				nInvalid++;
				bov_order_partial_update(invalid_order, tabInvalid, 0, nInvalid, 0);
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
		if (DRAWING) {
			bov_particles_draw(window, unchecked_points, 0, nPoints);
			bov_points_draw_with_order(window, invalid_points, invalid_order, 0, nInvalid);
			bov_points_draw_with_order(window, valid_points, valid_order, 0, nValid);
			if (use_verlet)
				bov_points_draw_with_order(window, verlet_points, verlet_order, 0, nVerlet);
			if (USE_IMPROVED_METHOD)
				bov_points_draw_with_order(window, finished_points, finished_order, 0, i- (int)(this_thread_number * ((double)nPoints) / NTH));
			bov_line_loop_draw(window, pointset, 0, BOV_TILL_END);
			if (use_verlet)
				bov_lines_draw(window, gridPoint, 0, BOV_TILL_END);
			bov_points_draw(window, this_point, 0, nPoints);
			bov_points_draw(window, checking_point, 0, 1);
			if (use_verlet)
				bov_points_draw(window, verlet_circle_point, 0, nPoints);
			bov_points_draw(window, circle_point, 0, nPoints);
			bov_window_update(window);
		}
		j_check = j;
		if (USE_IMPROVED_METHOD && f_check != frameCount / div) {
			j++;
		}
		else if (!USE_IMPROVED_METHOD)
			j = (frameCount / div) % (nPoints - 1) + (int)(i <= ((frameCount / div) % (nPoints - 1)));
		f_check = frameCount / div;
		frameCount++;
	}

	if (DRAWING) {
		free(tabValid);
		free(tabInvalid);
		free(tabFinished);
		bov_order_delete(valid_order);
		bov_order_delete(invalid_order);
		bov_order_delete(finished_order);
		bov_points_delete(this_point);
		bov_points_delete(circle_point);
		bov_points_delete(checking_point);
		bov_points_delete(valid_points);
		bov_points_delete(invalid_points);
		bov_points_delete(finished_points);
		bov_points_delete(unchecked_points);
		if (use_cells)
			bov_points_delete(gridPoint);
	}
	protected_inc(finished_threads);
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

//Fills the data vector
static void fillData(GLfloat(*data)[8], GLfloat(*coord)[2])
{
	float rmax = 100.0*sqrtf(2.0f);
	for(int i=0; i<NPTS; i++) {
		data[i][0] = rand() * 200.0 / RAND_MAX - 100.0; // x (rand between -100 and 100)
		data[i][1] = rand() * 200.0 / RAND_MAX - 100.0; // y (rand between -100 and 100)
		coord[i][0] = data[i][0];
		coord[i][1] = data[i][1];
		float r = sqrt(data[i][0] * data[i][0] + data[i][1] * data[i][1]);
		//data[i][2] = 0; //speed x (not used by default visualization)
		//data[i][3] = 0; // speed y (not used by default visualization)
		data[i][2] = rand() * 2.0 / RAND_MAX - 1.0; //Random starting speed
		data[i][3] = rand() * 2.0 / RAND_MAX - 1.0; //Random starting speed
		colormap(r/rmax, &data[i][4]); // fill color
		data[i][7] = 0.8f; // transparency
	}
}

//Changes the particle velocities randomly and updates the positions based, we assume elastic collisions with boundaries:
//	-timestep: time intervals at which these are updated
//	-xmin,xmax,ymin,ymax: boundaries of the domain
//	-maxspeed: the maximum speed that can be reached by the particles

void bouncyrandomupdate(GLfloat data[][8], GLfloat coord[][2],double timestep, double xmin, double xmax, double ymin, double ymax, double maxspeed) {
	for (int i = 0; i < NPTS; i++) {
		double speed = sqrtf(data[i][2] * data[i][2] + data[i][3] * data[i][3]);
		data[i][0] += data[i][2] * timestep;
		data[i][1] += data[i][3] * timestep;
		coord[i][0] = data[i][0];
		coord[i][1] = data[i][1];
		data[i][2] += ((double)rand() / RAND_MAX - 0.5)*(0.05*maxspeed) * timestep;
		data[i][3] += ((double)rand() / RAND_MAX - 0.5)*(0.05*maxspeed) * timestep;

		if (speed > maxspeed) {//Slows down if speed too high
			data[i][2] = data[i][2] * 0.9;
			data[i][3] = data[i][3] * 0.9;
		}

		//This next part of the code handles the cases where a particle bounces of the wall

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
	create_neighborhood_with_drawing(MAX_ITER);

	free(data);

	return EXIT_SUCCESS;
}