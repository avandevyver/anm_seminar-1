#ifndef NEIGHBORHOOD_SEARCH_H
#define NEIGHBORHOOD_SEARCH_H

#include "BOV.h"
#include <time.h>
#include <math.h>

extern int NPTS;

// see stringification process
#define xstr(s) str(s)
#define str(s) #s
#define M_PI 3.14159265358979323846

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
// index : the index in the data table, supposed to be available everywhere it is needed
// nNeighbours : number of neighbours in the linked list list
// nPotentialNeighbours : number of neighbours in the linked list potential_list, which contains the points to be checked in the verlet algorithm
// list : linked list of the actual neighbours of the particle at the index index in the data table
// potential_list : only list to be checked when we look for neighbours in the verlet algorithm
typedef struct neighborhood {
	int index;
	int nNeighbours;
	int nPotentialNeighbours;
	neighbours* list;
	neighbours* potential_list;
}neighborhood;

// Structure to be passed as argument to the function loop_without_drawing, now basically the same as the loop_arg_with_drawing without some useless parameters
// cells : array of size (size*size) that contains the cells of type cell
// cellCounter : counter to inform how many cells are and have been read already
// nPoints : number of points in the simulation
// iterations : used in the verlet algorithm; is equal to 0 each time we need to update the potential_list of the neighborhoods
// size : number of cells in a row
// nh : list of the neighborhoods to be filled; there are nPoints neighborhoods
// kh : size of the radius of the influence circle of a particle
// L : distance to be added to kh in the verlet algorithm; potential neighbours are the ones inside of a circle of radius kh+L
// coord : matrix of nPoints row and 2 column representing the positions of the particles used to draw; coord[i][0] == data[i][0] && coord[i][1] == data[i][1]
// data : matrix of nPoints row and 8 column representing the positions, speed, color and transparency of a particle
// use_verlet : int used as a boolean to inform if the verlet algorithm is used or not
// use_cells : int used as a boolean to inform if the cells are used or not
// use_improved_method : int used as a boolean to inform if the improved algorithm is used or not
typedef struct neighborhood_options {
	double kh;
	double L;
	int use_verlet;
	int use_cells;
	int use_improved_method;
	int half_length;
	int optimal_verlet_steps;
	neighborhood* nh;
}neighborhood_options;

// function used to print neighborhoods
// nh : array of neighborhoods to be printed
// data : table of the data's of the particles
// size : number of particles in the simulation
void printNeighborhood(neighborhood* nh, GLfloat(* data)[8]);

// function that basically fills the neighborhoods of the particles of one iteration, with the arguments args of type loop_arg
void neighborhood_update(neighborhood_options* options, neighborhood* nh, GLfloat(* data)[8], int iterations);

// function to change the particles velocities randomly and updates the positions based, we assume elastic collisions with boundaries
// data : table that contains the informations of the particles of the simulation
// coord : table that contains the positions of the particles of the simulations used to draw
// timestep : time intervals at which these are updated
// xmin,xmax,ymin,ymax : boundaries of the domain
// maxspeed : the maximum speed that can be reached by the particles
void bouncyrandomupdate(GLfloat(* data)[8], double timestep, double half_length, double maxspeed);

neighborhood_options* neighborhood_options_init(double timestep, double maxspeed);

void neighborhood_options_delete(neighborhood_options* options, neighborhood* nh);

int compare_neighborhoods(neighborhood* nh_1, neighborhood* nh_2);

#endif
