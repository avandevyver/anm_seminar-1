#include "neighborhood_search.h"
#include <math.h>


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

// function to create a neighbours to the neighborhood neigh[i] owned by the particles represented in data[i]
// index : the index of the neighbours to be added
// neigh : array of the neighborhood of the current iteration
// i : index if the particles that owns the neighborhood to modify
// d : distance between the particle in data[i] and the particle in data[index]
// is_after : ensures that all the neighbours added in the potential_list are after the owner of the neighborhood to modify; necessary when using improved algorithm
// is_potential : used to differentiate an actual neighbour and an only potential neighbour
void neighbours_new(int index, neighborhood* neigh, int i, double d, int is_after, int is_potential)
{
	if (is_after) {
		neighbours* newP = malloc(sizeof(neighbours));
		CHECK_MALLOC(newP);
		newP->index = index;
		newP->distance = d;
		newP->next = neigh[i].potential_list;
		neigh[i].potential_list = newP;
		neigh[i].nPotentialNeighbours++;
	}
	if (!is_potential) {
		neighbours* new = malloc(sizeof(neighbours));
		CHECK_MALLOC(new);
		new->index = index;
		new->distance = d;
		new->next = neigh[i].list;
		neigh[i].list = new;
		neigh[i].nNeighbours++;
	}
}

// function that frees the memory of a neighbour n passed as argument
void neighbours_delete(neighbours* n) {
	if (n) {
		neighbours* temp = n->next;
		free(n);
		neighbours_delete(temp);
	}
}


// function to create neighborhoods
// n : number of points, and so number of neighborhoods to create
// previous : pointer to the neighborhoods of the previous iteration, used in the verlet algorithm to keep the same potential_list when this should not be updated
neighborhood* neighborhood_new(neighborhood* previous, int iterations)
{
	for (int i = 0; i < NPTS; i++) {
			previous[i].nNeighbours = 0;
			neighbours_delete(previous[i].list);
			previous[i].list = NULL;
			if (!iterations) {
				previous[i].index = i;
				previous[i].nPotentialNeighbours = 0;
				neighbours_delete(previous[i].potential_list);
				previous[i].potential_list = NULL;
			}
	}
	return previous;
}

// function to properly delete the array of neighborhood nh, of size n
void neighborhood_delete(neighborhood* nh) {
	if (nh) {
		for (int i = 0; i < NPTS; i++) {
			neighbours_delete(nh[i].list);
			neighbours_delete(nh[i].potential_list);
		}
		free(nh);
	}
}

// function to create a node
// c : arrays of the cells of the simulation
// i : c[i] is the cell where the node should be included
// index : index of the particle to be added in the cell c[i]
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

// function to properly delete the node n given as argument
void node_delete(node* n) {
	if (n) {
		neighbours* temp = n->next;
		free(n);
		node_delete(temp);
	}
}

// function to create an array of cells of size n
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

// function to properly delete the array of cells c of size n given as arguments
void cell_delete(cell* c, GLsizei n) {
	if (c) {
		for (int i = 0; i < n; i++)
			node_delete(c[i].ResidentList);
		free(c);
	}
}

void printNeighborhood(neighborhood* nh, GLfloat(* data)[8]) {
	for (int i = 0; i < NPTS; i++) {
		printf("Resident %i : coordinate: %f %f   number of neighbours %i\n", i + 1, data[i][0], data[i][1], nh[i].nNeighbours);
		int j = 1;
		for (neighbours* current = nh[i].list; current; current = current->next)
			printf("   Neighbours %i : %f %f\n", j++, data[current->index][0], data[current->index][1]);
	}
}


// function used to print cells
// data : table of the data's of the particles
// c : array of cells to be printed
// size : number of cells in the simulation
void printCell(GLfloat(* data)[8], cell* c, int size) {
	for (int i = 0; i < size; i++) {
		printf("Cell %i : %i\n", i + 1, c[i].nResident);
		int j = 1;
		for (node* current = c[i].ResidentList; current; current = current->next)
			printf("   Neighbours %i : %f %f\n", j++, data[current->index][0], data[current->index][1]);
	}
}


// funtion to compute the radius kh of the circle of influence of a particle
// nPoints : number of particles in the simulation
// RA : int used as a boolean to choose over the algorithm of the radius choice; 
//      0 means the dummy algorithm, 1 means the more sophisticated algorithm expalined at the seminar ( (intersection between areaGrid and 1/4 areaCircle)/areaGrid = 21/nPoints)
double compute_kh(int RA) {
	double target = 21.0 / NPTS;
	if (NPTS < 21)
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


void neighborhood_update(neighborhood_options* options, neighborhood* nh, GLfloat(* data)[8], int iterations) {
	if (options->use_verlet)
		iterations = iterations % options->optimal_verlet_steps;
	else
		iterations = 0;

	if (options->use_verlet && iterations)
		neighborhood_new(nh,1);
	else
		neighborhood_new(nh,0);

	double kh = options->kh;
	int use_verlet = options->use_verlet;
	double L = 0.0;
	if (use_verlet) {
		L = options->L;
	}
	int use_improved_method = options->use_improved_method;
	int half_length = options->half_length;
	int use_cells = options->use_cells && (kh + L) < half_length / 3.0;
	int size = ceil(half_length / (kh + L));
	cell* cellArray = NULL;
	if (use_cells) {
		cellArray = cell_new(ceil(size) * ceil(size));
		for (int i = 0; i < NPTS; i++) {
			int cellNumber = ((int)((data[i][1] + half_length) / (2 * half_length) * size) * ceil(size) + (int)((data[i][0] + half_length) / (2 * half_length) * size));
			if (data[i][1] == half_length)
				cellNumber -= size;
			if (data[i][0] == half_length)
				cellNumber -= 1;
			node_new(cellArray, cellNumber, i);
		}
	}
	int cellCounter = 0;
	unsigned long frameCount = 0;
	int i = 0;
	int	j = 1;
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
	int are_still_neighbours = 1;
	while ((((use_verlet && iterations && use_cells) || !use_cells) && i < NPTS) || (use_cells && this_cell_number < size * size)) {
		if (i != i_check) {
			int are_still_neighbours = 1;
			if (use_verlet && iterations) {
				if (nh[i].potential_list) {
					checking_neighbours = nh[i].potential_list[0];
				}
				else {
					are_still_neighbours = 0;
				}
			}
			else if (use_cells) {
				if (i == 0) {
					this_cell_number = cellCounter;
					cellCounter++;
					while (!cellArray[this_cell_number].nResident && this_cell_number < size * size) {
						this_cell_number = cellCounter;
						cellCounter++;
					}
					if (this_cell_number < size * size) {
						this_cell = cellArray[this_cell_number];
						this_node = this_cell.ResidentList[0];
					}
				}
				else {
					if (this_node.next)
						this_node = *this_node.next;
					else {
						this_cell_number = cellCounter;
						cellCounter++;
						while (this_cell_number < size * size && !cellArray[this_cell_number].nResident && !cellArray[this_cell_number].ResidentList) {
							this_cell_number = cellCounter;
							cellCounter++;
						}
						if (this_cell_number < size * size) {
							this_cell = cellArray[this_cell_number];
							this_node = this_cell.ResidentList[0];
						}
					}
				}
				checked_cells = 0;
				if (use_improved_method && this_cell_number < size * size) {
					if (this_node.next) {
						checking_cell_number = this_cell_number;
						checking_cell = cellArray[checking_cell_number];
						checking_node = *this_node.next;
					}
					else {
						checking_cell_number = find_next_cell(this_cell_number, ++checked_cells, size, use_improved_method);
						while (checking_cell_number != -1 && !cellArray[checking_cell_number].nResident && !cellArray[checking_cell_number].ResidentList)
							checking_cell_number = find_next_cell(this_cell_number, ++checked_cells, size, use_improved_method);
						if (checking_cell_number != -1) {
							checking_cell = cellArray[checking_cell_number];
							checking_node = checking_cell.ResidentList[0];
						}
					}
				}
				else if (this_cell_number < size * size) {
					checking_cell_number = find_next_cell(this_cell_number, checked_cells, size, use_improved_method);
					while (checking_cell_number != -1 && !cellArray[checking_cell_number].nResident && !cellArray[checking_cell_number].ResidentList)
						checking_cell_number = find_next_cell(this_cell_number, ++checked_cells, size, use_improved_method);
					if (checking_cell_number != -1) {
						checking_cell = cellArray[checking_cell_number];
						checking_node = checking_cell.ResidentList[0];
					}
				}
			}
			if (use_improved_method && !(use_verlet && iterations)) {
				j_check = i;
				j = i + 1;
			}
		}
		i_check = i;
		if ((!(use_verlet && iterations) && use_cells && checking_cell_number == -1) || (use_verlet && iterations && !are_still_neighbours)) {
			i++;
			j_check = j;
		}
		if (j != j_check) {
			int index_i, index_j;
			if (use_verlet && iterations) {
				index_j = checking_neighbours.index;
				index_i = i;
			}
			else if (use_cells) {
				index_j = checking_node.index;
				index_i = this_node.index;
			}
			else {
				index_j = j;
				index_i = i;
			}
			double distance = sqrt((pow((double)data[index_j][0] - (double)data[index_i][0], 2) + pow((double)data[index_j][1] - (double)data[index_i][1], 2)));
			if (distance <= kh && index_i != index_j) {
				neighbours_new(index_j, nh, index_i, distance, !iterations, 0);
				if (use_improved_method)
					neighbours_new(index_i, nh, index_j, distance, 0, 0);
			}
			else if (use_verlet && !iterations && distance <= (kh + L) && index_i != index_j) {
				neighbours_new(index_j, nh, index_i, distance, !iterations, 1);
				if (use_improved_method)
					neighbours_new(index_i, nh, index_j, distance, 0, 1);
			}
			if (use_verlet && iterations) {
				if (checking_neighbours.next) {
					checking_neighbours = *(checking_neighbours.next);
				}
				else
					i++;
			}
			else if (use_cells)
				if (checking_node.next)
					checking_node = *(checking_node.next);
				else {
					checking_cell_number = find_next_cell(this_cell_number, ++checked_cells, size, use_improved_method);
					while (checking_cell_number != -1 && !cellArray[checking_cell_number].nResident && !cellArray[checking_cell_number].ResidentList)
						checking_cell_number = find_next_cell(this_cell_number, ++checked_cells, size, use_improved_method);
					if (checking_cell_number != -1) {
						checking_cell = cellArray[checking_cell_number];
						checking_node = cellArray[checking_cell_number].ResidentList[0];
					}
					else
						i++;
				}
			else
				if (j == NPTS - 1 || (i == NPTS - 1 && (j == NPTS - 2 || use_improved_method))) {
					i++;
				}
		}
		j_check = j;
		if (use_improved_method) {
			j++;
		}
		else
			j = (frameCount) % (NPTS - 1) + (int)(i <= ((frameCount) % (NPTS - 1)));
		frameCount++;
	}
	if (use_cells)
		cell_delete(cellArray, ceil(size) * ceil(size));
}

// function that returns which cell should be checked by a particle situated in this_cell
// this_cell : number of the cell that contains the particle which the neighborhood is filled
// checked_cell : number of cells that has already be completely checked for the current particle; serves as a first offset
// size : number of cells in one row, so there are (size*size) cells
// improve : serves as a second offset to avoid checking previous cells when the improved algorithm is activated
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

//Changes the particle velocities randomly and updates the positions based, we assume elastic collisions with boundaries:
//	-timestep: time intervals at which these are updated
//	-xmin,xmax,ymin,ymax: boundaries of the domain
//	-maxspeed: the maximum speed that can be reached by the particles

void bouncyrandomupdate(GLfloat(* data)[8], double timestep, double half_length, double maxspeed) {
	for (int i = 0; i < NPTS; i++) {
		double speed = sqrtf(data[i][2] * data[i][2] + data[i][3] * data[i][3]);
		data[i][0] += data[i][2] * timestep;
		data[i][1] += data[i][3] * timestep;
		data[i][2] += ((double)rand() / RAND_MAX - 0.5) * (0.05 * maxspeed) * timestep;
		data[i][3] += ((double)rand() / RAND_MAX - 0.5) * (0.05 * maxspeed) * timestep;

		if (speed > maxspeed) {//Slows down if speed too high
			data[i][2] = data[i][2] * 0.9;
			data[i][3] = data[i][3] * 0.9;
		}

		//This next part of the code handles the cases where a particle bounces of the wall

		//Particle is too high in x
		if (data[i][0] >= half_length) {
			data[i][0] -= 2 * (data[i][0] - half_length);
			data[i][2] = -data[i][2];
		}
		//Particle is too low in x
		if (data[i][0] <= -half_length) {
			data[i][0] -= 2 * (data[i][0] + half_length);
			data[i][2] = -data[i][2];
		}
		//Particle is too high in y
		if (data[i][1] >= half_length) {
			data[i][1] -= 2 * (data[i][1] - half_length);
			data[i][3] = -data[i][3];
		}
		//Particle is too low in y
		if (data[i][1] <= -half_length) {
			data[i][1] -= 2 * (data[i][1] + half_length);
			data[i][3] = -data[i][3];
		}
	}
}

// function that boils down to solving a cubic function and to find the optimal number of iterations without any update of the potential_list of the neighborhoods
// timestep : time intervals at which these are updated
// maxspeed : the maximum speed that can be reached by the particles
// kh : size of the radius of influence of a particle
int compute_optimal_verlet(double timestep, double maxspeed, double kh) {
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
		else {
			if (ceil(max_ans > floor(max_ans)))
				return ceil(max_ans);
			else
				return floor(max_ans);
		}
	}
}

neighborhood_options* neighborhood_options_init(double timestep, double maxspeed){
	neighborhood_options* options = malloc(sizeof(neighborhood_options));
	CHECK_MALLOC(options);

	int radius_algorithm = 1;

	options->half_length = 100;
	options->use_cells = 1;
	options->use_improved_method = 1;
	options->use_verlet = 1;
	options->optimal_verlet_steps = 0;
	options->kh = compute_kh(radius_algorithm) * 2 * options->half_length;
	options->L = 0.0;
	options->optimal_verlet_steps = compute_optimal_verlet(timestep, maxspeed, options->kh);
	if (options->optimal_verlet_steps == -1)
		options->use_verlet = 0;
	else
		options->L = options->optimal_verlet_steps * timestep * maxspeed;
	options->nh = calloc(NPTS, sizeof(neighborhood));
	CHECK_MALLOC(options->nh);
	return options;
}

void neighborhood_options_delete(neighborhood_options* options, neighborhood* nh) {
	if (options && options->nh && nh != options->nh)
		neighborhood_delete(options->nh);
	if (nh)
		neighborhood_delete(nh);
	if(options)
		free(options);
}

// function to check the equality of the 2 arrays of neighborhoods nh_1 and nh_2 of size nPoints
int compare_neighborhoods(neighborhood* nh_1, neighborhood* nh_2) {
	for (int i = 0; i < NPTS; i++) {
		if (nh_1[i].nNeighbours != nh_2[i].nNeighbours)
			return 0;
		neighbours* current = NULL;
		for (current = nh_1[i].list; current; current = current->next) {
			int is_equal = 0;
			neighbours* temp = NULL;
			for (temp = nh_2[i].list; temp && !is_equal; temp = temp->next) {
				is_equal = current->index == temp->index;
			}
			if (!is_equal) {
				return 0;
			}
		}
	}
	return 1;
}
