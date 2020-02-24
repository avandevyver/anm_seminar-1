#include "BOV.h"
#include <time.h>
#include <math.h>

// see stringification process
#define xstr(s) str(s)
#define str(s) #s
#define M_PI 3.14159265358979323846

// parameters
// To use them, fill the function create_neighborhood or create_neighborhood_with_drawing as in the original main
#define NPTS 50
#define USE_CELLS 1
#define USE_IMPROVED_METHOD 0
#define RADIUS_ALGORITHM 1
#define SEARCH_NEIGHBORHOOD 1
#define DRAWING 1
#define USE_VERLET 1
#define USE_THREADS 1
#define NTHREADS 4
#define MAX_ITER 5
#define KH 0

// malloc verification
// To use after each call to malloc, calloc and realloc
#define CHECK_MALLOC(ptr) if((ptr)==NULL) { \
		BOV_ERROR_LOG(BOV_OUT_OF_MEM_ERROR, "Memory allocation failed"); \
		exit(EXIT_FAILURE); }

// Structure to represent a neighbours, which is a node of the neighborhoods linked lists
// index : the index in the data table, supposed to be available everywhere it is needed
// distance : distance between the particle at the index in data and the particle that owns the neighborhood
// next : pointer to the next neighbours of the linked list
typedef struct neighbours {
	int index;
	double distance;
	struct neighbours* next;
}neighbours;

// Structure to represent the neighborhood of a particle
// mutex : used to protect the variables nNeighbours and list when more than one thread is used
// index : the index in the data table, supposed to be available everywhere it is needed
// nNeighbours : number of neighbours in the linked list list
// nPotentialNeighbours : number of neighbours in the linked list potential_list, which contains the points to be checked in the verlet algorithm
// list : linked list of the actual neighbours of the particle at the index index in the data table
// potential_list : only list to be checked when we look for neighbours in the verlet algorithm
typedef struct neighborhood {
	pthread_mutex_t mutex;
	int index;
	int nNeighbours;
	int nPotentialNeighbours;
	neighbours* list;
	neighbours* potential_list;
}neighborhood;

// Structure to represent a node of the linked list ResidentList of a cell, used to represent the particles contained in a cell
// index : the index in the data table, supposed to be available everywhere it is needed
// next : pointer to the next particle of the linked list
typedef struct node {
	int index;
	struct node* next;
}node;

// Structure to represent a cell and the particles it contains
// nResident : number of particles contained in the linked list ResidentList
// ResidentList : linked list used to contain the index of the particles that are in the current cell
typedef struct cell {
	int nResident;
	node* ResidentList;
}cell;

// Structure to associate an int with a mutex to protect it when using threads, should be modified by using the protected functions
// var : the value of the protected int
// mutex : mutex associated to the int var to protect it
typedef struct protected_int {
	int var;
	pthread_mutex_t mutex;
}protected_int;

// Structure to be passed as argument to the function loop_thread_with_drawing
// cells : array of size (size*size) that contains the cells of type cell
// cellCounter : counter to inform the threads how many cells are and have been read already
// finished_thread : used to inform the main thread how many thread are done, without blocking it with the use of pthread_join, and being able to update the original window
// nPoints : number of points in the simulation
// scale : used to scale the particle field
// div : used to represent the duration between each iteration; an iteration will last approximatively div/50 second(s)
// iterations : used in the verlet algorithm; is equal to 0 each time we need to update the potential_list of the neighborhoods
// size : number of cells in a row
// nh : list of the neighborhoods to be filled; there are nPoints neighborhoods
// kh : size of the radius of the influence circle of a particle
// L : distance to be added to kh in the verlet algorithm; potential neighbours are the ones inside of a circle of radius kh+L
// coord : matrix of nPoints row and 2 column representing the positions of the particles used to draw; coord[i][0] == data[i][0] && coord[i][1] == data[i][1]
// data : matrix of nPoints row and 8 column representing the positions, speed, color and transparency of a particle
// this_thread_number : number of the thread in which loop_thread_with_drawing is executed
// nThreads : total number of computing threads in the simulation
// use_verlet : int used as a boolean to inform if the verlet algorithm is used or not
// use_cells : int used as a boolean to inform if the cells are used or not
// use_improved_method : int used as a boolean to inform if the improved algorithm is used or not
typedef struct loop_thread_arg_wiht_drawing {
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
	int nThreads;
	int use_verlet;
	int use_cells;
	int use_improved_method;
}ltad;

// Structure to be passed as argument to the function loop_thread, now basically the same as the loop_thread_arg_with_drawing without some useless parameters
// cells : array of size (size*size) that contains the cells of type cell
// cellCounter : counter to inform the threads how many cells are and have been read already
// nPoints : number of points in the simulation
// iterations : used in the verlet algorithm; is equal to 0 each time we need to update the potential_list of the neighborhoods
// size : number of cells in a row
// nh : list of the neighborhoods to be filled; there are nPoints neighborhoods
// kh : size of the radius of the influence circle of a particle
// L : distance to be added to kh in the verlet algorithm; potential neighbours are the ones inside of a circle of radius kh+L
// coord : matrix of nPoints row and 2 column representing the positions of the particles used to draw; coord[i][0] == data[i][0] && coord[i][1] == data[i][1]
// data : matrix of nPoints row and 8 column representing the positions, speed, color and transparency of a particle
// this_thread_number : number of the thread in which loop_thread is executed
// nThreads : total number of computing threads in the simulation
// use_verlet : int used as a boolean to inform if the verlet algorithm is used or not
// use_cells : int used as a boolean to inform if the cells are used or not
// use_improved_method : int used as a boolean to inform if the improved algorithm is used or not
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
	int nThreads;
	int use_verlet;
	int use_cells;
	int use_improved_method;
}lta;

// function to create a neighbours to the neighborhood neigh[i] owned by the particles represented in data[i]
// index : the index of the neighbours to be added
// neigh : array of the neighborhood of the current iteration
// i : index if the particles that owns the neighborhood to modify
// d : distance between the particle in data[i] and the particle in data[index]
// is_after : ensures that all the neighbours added in the potential_list are after the owner of the neighborhood to modify; necessary when using improved algorithm
// is_potential : used to differentiate an actual neighbour and an only potential neighbour
void neighbours_new(int index, neighborhood* neigh, int i, float d, int is_after, int is_potential);

// function that frees the memory of a neighbour n passed as argument
void neighbours_delete(neighbours* n);

// function to create neighborhoods
// n : number of points, and so number of neighborhoods to create
// previous : pointer to the neighborhoods of the previous iteration, used in the verlet algorithm to keep the same potential_list when this should not be updated
neighborhood* neighborhood_new(GLsizei n, neighborhood* previous);

// function to properly delete the array of neighborhood nh, of size n
void neighborhood_delete(neighborhood* nh, GLsizei n);

// function to create a node
// c : arrays of the cells of the simulation
// i : c[i] is the cell where the node should be included
// index : index of the particle to be added in the cell c[i]
struct node* node_new(cell* c, int i, int index);

// function to properly delete the node n given as argument
void node_delete(node* n);

// function to create an array of cells of size n
cell* cell_new(GLsizei n);

// function to properly delete the array of cells c of size n given as arguments
void cell_delete(cell* c, GLsizei n);

// function used to print neighborhoods
// nh : array of neighborhoods to be printed
// data : table of the data's of the particles
// size : number of particles in the simulation
void printNeighborhood(neighborhood* nh, GLfloat data[][8], int size);

// function used to print cells
// data : table of the data's of the particles
// c : array of cells to be printed
// size : number of cells in the simulation
void printCell(GLfloat data[][8], cell* c, int size);

// for v, a floating point value between 0 and 1, this function fills color with
// the improved jet colormap color corresponding to v
static void colormap(float v, float color[3]);

// function used to safely increment the value var of the protected int guard given as argument
void protected_inc(protected_int* guard);

// function used to safely get the value var of the protected int guard given as argument
int protected_get(protected_int* guard);

// function used to safely get and then increment the value var of the protected int guard given as argument
int protected_get_inc(protected_int* guard);

// function used to safely reset the value var (setting it to 0) of the protected int guard given as argument
void protected_reset(protected_int* guard);

// funtion to compute the radius kh of the circle of influence of a particle
// nPoints : number of particles in the simulation
// RA : int used as a boolean to choose over the algorithm of the radius choice; 
//      0 means the dummy algorithm, 1 means the more sophisticated algorithm expalined at the seminar ( (intersection between areaGrid and 1/4 areaCircle)/areaGrid = 21/nPoints)
double compute_kh(int nPoints, int RA);

// function to find the points that will be bound to form the cells on the drawing
// kh : size of the radius of influence of a particle, and hence of the width and length of a cell
// xmax : half size of the grid's length
bov_points_t* find_grid_points(double kh, double xmax);

// function that basically fills the neighborhoods of the particles of one iteration, with the arguments args of type loop_thread_arg
void* loop_thread(void* args);

// function that basically fills the neighborhoods of the particles of one iteration, as the previous one
// but with the arguments args of type loop_thread_arg_with_drawing, and of course with the drawing of the main grid and each steps
void* loop_thread_with_drawing(void* args);

// function that returns which cell should be checked by a particle situated in this_cell
// this_cell : number of the cell that contains the particle which the neighborhood is filled
// checked_cell : number of cells that has already be completely checked for the current particle; serves as a first offset
// size : number of cells in one row, so there are (size*size) cells
// improve : serves as a second offset to avoid checking previous cells when the improved algorithm is activated
int find_next_cell(int this_cell, int checked_cells, int size, int improve);

// function to fill the data table of the nPoints particles positions, speeds, colors and transparency and the coord table with the nPoints particles positions used to draw;
// data[i][0] == coord[i][0] && data[i][1] == coord[i][1]
void fillData(GLfloat(*data)[8], GLfloat(*coord)[2], int nPoints);

// function to change the particles velocities randomly and updates the positions based, we assume elastic collisions with boundaries
// data : table that contains the informations of the particles of the simulation
// coord : table that contains the positions of the particles of the simulations used to draw
// timestep : time intervals at which these are updated
// xmin,xmax,ymin,ymax : boundaries of the domain
// maxspeed : the maximum speed that can be reached by the particles
void bouncyrandomupdate(GLfloat(*data)[8], GLfloat(*coord)[2], double timestep, double xmin, double xmax, double ymin, double ymax, double maxspeed);

// function that boils down to solving a cubic function and to find the optimal number of iterations without any update of the potential_list of the neighborhoods
// timestep : time intervals at which these are updated
// maxspeed : the maximum speed that can be reached by the particles
// kh : size of the radius of influence of a particle
int compute_optimal_verlet(double timestep, double maxspeed, double kh);

// main function of the neighborhood search algorithm that return the list of array of neighborhood of the iter number of iterations and draw the steps
// data : table that contains the informations of the particles of the simulation
// coord : table that contains the positions of the particles of the simulations used to draw
// iter : number of iterations of the complete neighborhoods search, particles are moved after each iterations and neighborhoods are then calculated
// data_filled : int used as a boolean to warn if the tables data and coord has to be filled or not
// nPoints : number of particles in the simulation
// nThreads : total number of computing threads in the simulation
// use_cells : int used as a boolean to inform if the cells are used or not
// use_improved_method : int used as a boolean to inform if the improved algorithm is used or not
// radius_algorithm : int used as a boolean to choose over the algorithm of the radius choice; 
//                    0 means the dummy algorithm, 1 means the more sophisticated algorithm expalined at the seminar ( (intersection between areaGrid and 1/4 areaCircle)/areaGrid = 21/nPoints)
// search_algorithm : int used as a boolean to choose if the neighborhood search has to be performed or not, allows a faster visualisation of the particles movement
// use_verlet : int used as a boolean to inform if the verlet algorithm is used or not
// use_threads :  int used as a boolean to inform if threads are used or not
// nThreads : total number of computing threads in the simulation
// seed : used to seed the random; with the same seed, the same data and coord are filled, and the same iterations occured; useful for example to compare 2 different algorithms results
neighborhood** create_neighborhood_with_drawing(GLfloat(*data)[8], GLfloat(*coord)[2], int iter, int data_filled, int nPoints, int use_cells, int use_improved_method, int radius_algorithm, int search_neighborhood, int use_verlet, int use_threads, int nThreads, unsigned int seed);

// main function of the neighborhood search algorithm that return the list of array of neighborhood of the iter number of iterations without drawing the steps
// data : table that contains the informations of the particles of the simulation
// coord : table that contains the positions of the particles of the simulations used to draw
// iter : number of iterations of the complete neighborhoods search, particles are moved after each iterations and neighborhoods are then calculated
// data_filled : int used as a boolean to warn if the tables data and coord has to be filled or not
// nPoints : number of particles in the simulation
// nThreads : total number of computing threads in the simulation
// use_cells : int used as a boolean to inform if the cells are used or not
// use_improved_method : int used as a boolean to inform if the improved algorithm is used or not
// radius_algorithm : int used as a boolean to choose over the algorithm of the radius choice; 
//                    0 means the dummy algorithm, 1 means the more sophisticated algorithm expalined at the seminar ( (intersection between areaGrid and 1/4 areaCircle)/areaGrid = 21/nPoints)
// search_algorithm : int used as a boolean to choose if the neighborhood search has to be performed or not, allows a faster visualisation of the particles movement
// use_verlet : int used as a boolean to inform if the verlet algorithm is used or not
// use_threads :  int used as a boolean to inform if threads are used or not
// nThreads : total number of computing threads in the simulation
// seed : used to seed the random; with the same seed, the same data and coord are filled, and the same iterations occured; useful for example to compare 2 different algorithms results
neighborhood** create_neighborhood(GLfloat(*data)[8], GLfloat(*coord)[2], int iter, int data_filled, int nPoints, int use_cells, int use_improved_method, int radius_algorithm, int use_verlet, int use_threads, int nThreads, unsigned int seed);

// function to check the equality of the 2 lists of arrays of neighborhoods nhList_1 and nhList_2;
// the size of the list is given by iter and the size of the neighborhoods is given by nPoints
int compare_neighborhood_lists(neighborhood** nhList_1, neighborhood** nhList_2, int iter, int nPoints);

// function to check the equality of the 2 arrays of neighborhoods nh_1 and nh_2 of size nPoints
int compare_neighborhoods(neighborhood* nh_1, neighborhood* nh_2, int nPoints);

// function only used to show the equivalence of some of the major algorithm combination, with our without drawing the steps
void show_equivalence_of_major_algorithm_combination(int drawing);