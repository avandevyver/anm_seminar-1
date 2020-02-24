#include "neighborhood_search_for_mac.h"

void neighbours_new_for_mac(int index, neighborhood_for_mac* neigh, int i, double d, int is_after, int is_potential)
{
	if (is_after) {
		neighbours_for_mac* newP = malloc(sizeof(neighbours_for_mac));
		CHECK_MALLOC(newP);
		newP->index = index;
		newP->distance = d;
		newP->next = neigh[i].potential_list;
		neigh[i].potential_list = newP;
		neigh[i].nPotentialNeighbours++;
	}
	if (!is_potential) {
		neighbours_for_mac* new = malloc(sizeof(neighbours_for_mac));
		CHECK_MALLOC(new);
		new->index = index;
		new->distance = d;
		new->next = neigh[i].list;
		neigh[i].list = new;
		neigh[i].nNeighbours++;
	}
}

void neighbours_delete_for_mac(neighbours_for_mac* n) {
	if (n) {
		neighbours_for_mac* temp = n->next;
		free(n);
		neighbours_delete_for_mac(temp);
	}
}

neighborhood_for_mac* neighborhood_new_for_mac(GLsizei n, neighborhood_for_mac* previous)
{
	neighborhood_for_mac* neigh = calloc(n, sizeof(neighborhood_for_mac));
	CHECK_MALLOC(neigh);
	for (int i = 0; i < n; i++) {
		neigh[i].index = i;
		neigh[i].nNeighbours = 0;
		neigh[i].list = NULL;
		if (previous) {
			neigh[i].nPotentialNeighbours = previous[i].nPotentialNeighbours;
			neigh[i].potential_list = previous[i].potential_list;
		}
		else {
			neigh[i].nPotentialNeighbours = 0;
			neigh[i].potential_list = NULL;
		}
	}
	return neigh;
}

void neighborhood_delete_for_mac(neighborhood_for_mac* nh, GLsizei n) {
	for (int i = 0; i < n; i++) {
		neighbours_delete_for_mac(nh[i].list);
		//neighbours_delete(nh[i].potential_list);
	}
	free(nh);
}

struct node_for_mac* node_new_for_mac(cell_for_mac* c, int i, int index)
{
	node_for_mac* new = malloc(sizeof(node_for_mac));
	CHECK_MALLOC(new);
	new->index = index;
	new->next = c[i].ResidentList;
	c[i].ResidentList = new;
	c[i].nResident++;
	return new;
}

void node_delete_for_mac(node_for_mac* n) {
	if (n) {
		neighbours_for_mac* temp = n->next;
		free(n);
		node_delete_for_mac(temp);
	}
}

cell_for_mac* cell_new_for_mac(GLsizei n)
{
	cell_for_mac* c = calloc(n, sizeof(cell_for_mac));
	CHECK_MALLOC(c);
	for (int i = 0; i < n; i++) {
		c[i].nResident = 0;
		c[i].ResidentList = NULL;
	}
	return c;
}

void cell_delete_for_mac(cell_for_mac* c, GLsizei n) {
	for (int i = 0; i < n; i++)
		node_delete_for_mac(c[i].ResidentList);
	free(c);
}

void printNeighborhood_for_mac(neighborhood_for_mac* nh, GLfloat data[][8], int size) {
	for (int i = 0; i < size; i++) {
		printf("Resident %i : coordinate: %f %f   number of neighbours %i\n", i + 1, data[i][0], data[i][1], nh[i].nNeighbours);
		int j = 1;
		for (neighbours_for_mac* current = nh[i].list; current; current = current->next)
			printf("   Neighbours %i : %f %f\n", j++, data[current->index][0], data[current->index][1]);
	}
}

void printCell_for_mac(GLfloat data[][8], cell_for_mac* c, int size) {
	for (int i = 0; i < size; i++) {
		printf("Cell %i : %i\n", i + 1, c[i].nResident);
		int j = 1;
		for (node_for_mac* current = c[i].ResidentList; current; current = current->next)
			printf("   Neighbours %i : %f %f\n", j++, data[current->index][0], data[current->index][1]);
	}
}

// for v, a floating point value between 0 and 1, this function fills color with
// the improved jet colormap color corresponding to v
static void colormap_for_mac(float v, float color[3])
{
	float v1 = 3.5 * (v - 0.7);
	float v2 = 1.25 * v;
	float v3 = fminf(0.5, v) * 2.0;

	color[0] = -v1 * v1 + 1.0f;
	color[1] = 6.0f * v2 * v2 * (1.0f - v2);
	color[2] = 5.5f * v3 * (1.0f - v3) * (1.0f - v3);

	// alternative: classical jet colormap
	// color[0] = 1.5 - 4.0 * fabs(v - 0.75);
	// color[1] = 1.5 - 4.0 * fabs(v - 0.5 );
	// color[2] = 1.5 - 4.0 * fabs(v - 0.25);
}

double compute_kh_for_mac(int nPoints, int RA) {
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

bov_points_t* find_grid_points_for_mac(double kh, double xmax) {
	double size = 2 * xmax / kh;
	float(*pts)[2] = malloc((ceil(size) - 1) * 4 * sizeof(pts[0]));
	CHECK_MALLOC(pts);
	for (int i = 1; i < size; i++) {
		pts[(i - 1) * 4][0] = xmax * (2.0 * i / (size)-1.0);
		pts[(i - 1) * 4][1] = xmax * (-1.0);
		pts[(i - 1) * 4 + 1][0] = xmax * (2.0 * i / (size)-1.0);
		pts[(i - 1) * 4 + 1][1] = xmax * (1.0);
		pts[(i - 1) * 4 + 2][0] = xmax * (-1.0);
		pts[(i - 1) * 4 + 2][1] = xmax * (2.0 * i / (size)-1.0);
		pts[(i - 1) * 4 + 3][0] = xmax * (1.0);
		pts[(i - 1) * 4 + 3][1] = xmax * (2.0 * i / (size)-1.0);
	}
	bov_points_t* points = bov_points_new(pts, (ceil(size) - 1) * 4, GL_STATIC_DRAW);
	return points;
}


void* loop_without_drawing_for_mac(void* args) {
	lafm* loop = (lafm*)args;
	cell_for_mac* cellArray = loop->cells;
	int cellCounter = loop->cellCounter;
	int size = ceil(loop->size);
	neighborhood_for_mac* nh = loop->nh;
	double kh = loop->kh;
	double L = loop->L;
	GLsizei nPoints = loop->nPoints;
	int iterations = loop->iterations;
	int use_verlet = loop->use_verlet;
	int use_cells = loop->use_cells;
	int use_improved_method = loop->use_improved_method;
	unsigned long frameCount = 0;
	int i, j;
	i = 0;
	if (use_improved_method)
		j = i + 1;
	else
		j = 0;
	int nValid = 0;
	int nInvalid = 0;
	int nVerlet = 0;
	int checking_cell_number = -1;
	int this_cell_number = -1;
	int checked_cells = 0;
	cell_for_mac checking_cell;
	node_for_mac checking_node;
	neighbours_for_mac checking_neighbours;
	cell_for_mac this_cell;
	node_for_mac this_node;
	int i_check = i - 1;
	int j_check = j - 1;
	int f_check = 1;
	int are_still_neighbours = 1;
	while (use_cells && !(use_verlet && iterations) && i < nPoints && this_cell_number < size * size) {
		if (i != i_check) {
			nValid = 0;
			nInvalid = 0;
			nVerlet = 0;
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
						checking_cell_number = find_next_cell_for_mac(this_cell_number, ++checked_cells, size, use_improved_method);
						while (checking_cell_number != -1 && !cellArray[checking_cell_number].nResident && !cellArray[checking_cell_number].ResidentList)
							checking_cell_number = find_next_cell_for_mac(this_cell_number, ++checked_cells, size, use_improved_method);
						if (checking_cell_number != -1) {
							checking_cell = cellArray[checking_cell_number];
							checking_node = checking_cell.ResidentList[0];
						}
					}
				}
				else if (this_cell_number < size * size) {
					checking_cell_number = find_next_cell_for_mac(this_cell_number, checked_cells, size, use_improved_method);
					while (checking_cell_number != -1 && !cellArray[checking_cell_number].nResident && !cellArray[checking_cell_number].ResidentList)
						checking_cell_number = find_next_cell_for_mac(this_cell_number, ++checked_cells, size, use_improved_method);
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
			double distance = sqrt((pow((double)loop->coord[index_j][0] - (double)loop->coord[index_i][0], 2) + pow((double)loop->coord[index_j][1] - (double)loop->coord[index_i][1], 2)));
			if (distance <= kh && index_i != index_j) {
				neighbours_new_for_mac(index_j, nh, index_i, distance, !iterations, 0);
				if (use_improved_method)
					neighbours_new_for_mac(index_i, nh, index_j, distance, 0, 0);
			}
			else if (use_verlet && !iterations && distance <= (kh + L) && index_i != index_j) {
				neighbours_new_for_mac(index_j, nh, index_i, distance, !iterations, 1);
				if (use_improved_method)
					neighbours_new_for_mac(index_i, nh, index_j, distance, 0, 1);
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
					checking_cell_number = find_next_cell_for_mac(this_cell_number, ++checked_cells, size, use_improved_method);
					while (checking_cell_number != -1 && !cellArray[checking_cell_number].nResident && !cellArray[checking_cell_number].ResidentList)
						checking_cell_number = find_next_cell_for_mac(this_cell_number, ++checked_cells, size, use_improved_method);
					if (checking_cell_number != -1) {
						checking_cell = cellArray[checking_cell_number];
						checking_node = cellArray[checking_cell_number].ResidentList[0];
					}
					else
						i++;
				}
			else
				if (j == nPoints - 1 || (i == nPoints - 1 && (j == nPoints - 2 || use_improved_method))) {
					i++;
				}
		}
		j_check = j;
		if (use_improved_method) {
			j++;
		}
		else if (!use_improved_method)
			j = (frameCount) % (nPoints - 1) + (int)(i <= ((frameCount) % (nPoints - 1)));
		frameCount++;
	}
}

void* loop_with_drawing_for_mac(void* args) {
	ladfm* loop = (ladfm*)args;
	cell_for_mac* cellArray = loop->cells;
	int cellCounter = loop->cellCounter;
	int div = loop->div;
	int size = ceil(loop->size);
	neighborhood_for_mac* nh = loop->nh;
	double kh = loop->kh;
	double L = loop->L;
	double scale = loop->scale;
	GLsizei nPoints = loop->nPoints;
	int iterations = loop->iterations;
	int use_verlet = loop->use_verlet;
	int use_cells = loop->use_cells;
	int use_improved_method = loop->use_improved_method;
	double points_width = 0.02;
	char name[32];
	snprintf(name, sizeof(name), "Neighborhood: loop");
	bov_window_t* window = NULL;
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
		CHECK_MALLOC(tabValid);
		tabInvalid = calloc(nPoints, sizeof(GLuint));
		CHECK_MALLOC(tabInvalid);
		if (use_improved_method) {
			tabFinished = calloc(nPoints, sizeof(GLuint));
			CHECK_MALLOC(tabFinished);
			finished_order = bov_order_new(tabFinished, nPoints + 1, GL_DYNAMIC_DRAW);
		}
		valid_order = bov_order_new(tabValid, nPoints + 1, GL_DYNAMIC_DRAW);
		invalid_order = bov_order_new(tabInvalid, nPoints + 1, GL_DYNAMIC_DRAW);
		unchecked_points = bov_particles_new(loop->data, nPoints, GL_DYNAMIC_DRAW);
		this_point = bov_points_new(loop->coord, 1, GL_DYNAMIC_DRAW);
		if (use_verlet) {
			tabVerlet = calloc(nPoints, sizeof(GLuint));
			CHECK_MALLOC(tabVerlet);
			verlet_circle_point = bov_points_new(loop->coord, 1, GL_DYNAMIC_DRAW);
			verlet_order = bov_order_new(tabVerlet, nPoints + 1, GL_DYNAMIC_DRAW);
			verlet_points = bov_points_new(loop->coord, nPoints, GL_STATIC_DRAW);
		}
		circle_point = bov_points_new(loop->coord, 1, GL_DYNAMIC_DRAW);
		checking_point = bov_points_new(loop->coord, 0, GL_DYNAMIC_DRAW);
		valid_points = bov_points_new(loop->coord, nPoints, GL_STATIC_DRAW);
		invalid_points = bov_points_new(loop->coord, nPoints, GL_STATIC_DRAW);
		if (use_improved_method)
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
			gridPoint = find_grid_points(kh + L, 100);
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
		bov_points_set_color(this_point, test_point_color);
		bov_points_set_color(checking_point, other_color);
		bov_points_set_color(circle_point, circle_color);
		bov_points_set_color(valid_points, valid_color);
		bov_points_set_color(invalid_points, invalid_color);
		if (use_improved_method)
			bov_points_set_color(finished_points, finished_color);
		bov_points_set_outline_width(circle_point, 0.);
		bov_points_set_width(circle_point, kh * scale);
		bov_points_set_width(this_point, points_width);
		bov_points_set_width(checking_point, points_width);
		bov_points_set_width(unchecked_points, points_width);
		bov_points_set_width(valid_points, points_width);
		bov_points_set_width(invalid_points, points_width);
		if (use_improved_method)
			bov_points_set_width(finished_points, points_width);
		bov_points_set_outline_width(this_point, 0.0025);
		bov_points_set_outline_width(checking_point, 0.0025);
		bov_points_set_outline_width(unchecked_points, 0.0025);
		bov_points_set_outline_width(valid_points, 0.0025);
		bov_points_set_outline_width(invalid_points, 0.0025);
		if (use_improved_method)
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
		if (use_improved_method) {
			bov_points_scale(finished_points, (GLfloat[2]) { scale, scale });
			bov_points_set_pos(finished_points, (GLfloat[2]) { 0.0, -0.1 });
		}
		bov_points_scale(circle_point, (GLfloat[2]) { scale, scale });
		bov_points_set_pos(circle_point, (GLfloat[2]) { 0.0, -0.1 });
	}
	unsigned long frameCount = 0;
	int i, j;
	i = 0;
	if (use_improved_method)
		j = i + 1;
	else
		j = 0;
	int nValid = 0;
	int nInvalid = 0;
	int nVerlet = 0;
	int checking_cell_number = -1;
	int this_cell_number = -1;
	int checked_cells = 0;
	cell_for_mac checking_cell;
	node_for_mac checking_node;
	neighbours_for_mac checking_neighbours;
	cell_for_mac this_cell;
	node_for_mac this_node;
	int i_check = i - 1;
	int j_check = j - 1;
	int f_check = 1;
	int are_still_neighbours = 1;
	if (!DRAWING)
		div = 1;
	while ((!DRAWING || !bov_window_should_close(window)) && i<nPoints && this_cell_number < size * size) {
		if (i != i_check) {
			nValid = 0;
			nInvalid = 0;
			nVerlet = 0;
			are_still_neighbours = 1;
			if (use_verlet && iterations) {
				if (nh[i].potential_list) {
					checking_neighbours = nh[i].potential_list[0];
				}
				else {
					are_still_neighbours = 0;
				}
				if (DRAWING) {
					bov_points_update(this_point, loop->coord[i], 1);
					bov_points_update(circle_point, loop->coord[i], 1);
					if (use_verlet)
						bov_points_update(verlet_circle_point, loop->coord[i], 1);
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
						checking_cell_number = find_next_cell_for_mac(this_cell_number, ++checked_cells, size, use_improved_method);
						while (checking_cell_number != -1 && !cellArray[checking_cell_number].nResident && !cellArray[checking_cell_number].ResidentList)
							checking_cell_number = find_next_cell_for_mac(this_cell_number, ++checked_cells, size, use_improved_method);
						if (checking_cell_number != -1) {
							checking_cell = cellArray[checking_cell_number];
							checking_node = checking_cell.ResidentList[0];
						}
					}
				}
				else if (this_cell_number < size * size) {
					checking_cell_number = find_next_cell_for_mac(this_cell_number, checked_cells, size, use_improved_method);
					while (checking_cell_number != -1 && !cellArray[checking_cell_number].nResident && !cellArray[checking_cell_number].ResidentList)
						checking_cell_number = find_next_cell_for_mac(this_cell_number, ++checked_cells, size, use_improved_method);
					if (checking_cell_number != -1) {
						checking_cell = cellArray[checking_cell_number];
						checking_node = checking_cell.ResidentList[0];
					}
				}
				if (DRAWING && this_cell_number < size * size) {
					bov_points_update(this_point, loop->coord[this_node.index], 1);
					bov_points_update(circle_point, loop->coord[this_node.index], 1);
					if (use_verlet)
						bov_points_update(verlet_circle_point, loop->coord[this_node.index], 1);
				}
			}
			else if (DRAWING) {
				bov_points_update(this_point, loop->coord[i], 1);
				bov_points_update(circle_point, loop->coord[i], 1);
				if (use_verlet)
					bov_points_update(verlet_circle_point, loop->coord[i], 1);
			}
			if (use_improved_method) {
				if (DRAWING) {
					if (use_cells && !(use_verlet && iterations))
						tabFinished[i] = this_node.index;
					else
						tabFinished[i] = i;
					finished_order = bov_order_partial_update(finished_order, tabFinished, 0, i, 0);
				}
				if (!(use_verlet && iterations)) {
					j_check = i;
					j = i + 1;
				}
			}
		}
		i_check = i;
		if ((use_cells && checking_cell_number == -1 && !(use_verlet && iterations)) || (use_verlet && iterations && !are_still_neighbours)) {
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
			if (DRAWING)
				bov_points_update(checking_point, loop->coord[index_j], 1);
			double distance = (pow((double)loop->coord[index_j][0] - (double)loop->coord[index_i][0], 2) + pow((double)loop->coord[index_j][1] - (double)loop->coord[index_i][1], 2));
			if (distance <= pow(kh, 2) && index_i != index_j) {
				if (DRAWING) {
					tabValid[nValid] = index_j;
					nValid++;
					bov_order_partial_update(valid_order, tabValid, 0, nValid, 0);
				}
				neighbours_new_for_mac(index_j, nh, index_i, distance, !iterations, 0);
				if (use_improved_method)
					neighbours_new_for_mac(index_i, nh, index_j, distance, 0, 0);
			}
			else if (use_verlet && !iterations && distance <= (kh + L) && index_i != index_j) {
				if (DRAWING) {
					tabVerlet[nVerlet] = index_j;
					nVerlet++;
					bov_order_partial_update(verlet_order, tabVerlet, 0, nVerlet, 0);
				}
				neighbours_new_for_mac(index_j, nh, index_i, distance, !iterations, 1);
				if (use_improved_method)
					neighbours_new_for_mac(index_i, nh, index_j, distance, 0, 1);
			}
			else if (DRAWING) {
				tabInvalid[nInvalid] = index_j;
				nInvalid++;
				bov_order_partial_update(invalid_order, tabInvalid, 0, nInvalid, 0);
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
					checking_cell_number = find_next_cell_for_mac(this_cell_number, ++checked_cells, size, use_improved_method);
					while (checking_cell_number != -1 && !cellArray[checking_cell_number].nResident && !cellArray[checking_cell_number].ResidentList)
						checking_cell_number = find_next_cell_for_mac(this_cell_number, ++checked_cells, size, use_improved_method);
					if (checking_cell_number != -1) {
						checking_cell = cellArray[checking_cell_number];
						checking_node = cellArray[checking_cell_number].ResidentList[0];
					}
					else
						i++;
				}
			else
				if (j == nPoints - 1 || (i == nPoints - 1 && (j == nPoints - 2 || use_improved_method))) {
					i++;
				}
		}
		if (DRAWING) {
			bov_particles_draw(window, unchecked_points, 0, BOV_TILL_END);
			bov_points_draw_with_order(window, invalid_points, invalid_order, 0, nInvalid);
			bov_points_draw_with_order(window, valid_points, valid_order, 0, nValid);
			if (use_verlet)
				bov_points_draw_with_order(window, verlet_points, verlet_order, 0, nVerlet);
			if (use_improved_method)
					bov_points_draw_with_order(window, finished_points, finished_order, 0, i);
			bov_line_loop_draw(window, pointset, 0, BOV_TILL_END);
			if (use_cells)
				bov_lines_draw(window, gridPoint, 0, BOV_TILL_END);
			bov_points_draw(window, this_point, 0, nPoints);
			bov_points_draw(window, checking_point, 0, 1);
			if (use_verlet)
				bov_points_draw(window, verlet_circle_point, 0, nPoints);
			bov_points_draw(window, circle_point, 0, nPoints);
			bov_window_update(window);
		}
		j_check = j;
		if (use_improved_method && f_check != frameCount / div) {
			j++;
		}
		else if (!use_improved_method)
			j = (frameCount / div) % (nPoints - 1) + (int)(i <= ((frameCount / div) % (nPoints - 1)));
		f_check = frameCount / div;
		frameCount++;
	}

	free(tabValid);
	free(tabInvalid);
	if (use_improved_method)
		free(tabFinished);
	bov_order_delete(valid_order);
	bov_order_delete(invalid_order);
	if (use_improved_method)
		bov_order_delete(finished_order);
	bov_points_delete(this_point);
	bov_points_delete(circle_point);
	bov_points_delete(checking_point);
	bov_points_delete(valid_points);
	bov_points_delete(invalid_points);
	if (use_improved_method)
		bov_points_delete(finished_points);
	bov_points_delete(unchecked_points);
	if (use_cells)
		bov_points_delete(gridPoint);
	return(window);
}

int find_next_cell_for_mac(int this_cell, int checked_cells, int size, int improve) {
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
void fillData_for_mac(GLfloat(*data)[8], GLfloat(*coord)[2], int nPoints)
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
		colormap_for_mac(r / rmax, &data[i][4]); // fill color
		data[i][7] = 0.8f; // transparency
	}
}

//Changes the particle velocities randomly and updates the positions based, we assume elastic collisions with boundaries:
//	-timestep: time intervals at which these are updated
//	-xmin,xmax,ymin,ymax: boundaries of the domain
//	-maxspeed: the maximum speed that can be reached by the particles

void bouncyrandomupdate_for_mac(GLfloat(*data)[8], GLfloat(*coord)[2], double timestep, double xmin, double xmax, double ymin, double ymax, double maxspeed) {
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

int compute_optimal_verlet_for_mac(double timestep, double maxspeed, double kh) {
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

neighborhood_for_mac** create_neighborhood_with_drawing_for_mac(GLfloat(*data)[8], GLfloat(*coord)[2], int iter, int data_filled, int nPoints, int use_cells, int use_improved_method, int radius_algorithm, int search_neighborhood, int use_verlet, unsigned int seed) {
	//Defining the domain
	int xmax = 100;
	int xmin = -100;
	int ymax = 100;
	int ymin = -100;
	bov_window_t* window = bov_window_new(800, 800, "ANM Project: SPH");
	bov_window_set_color(window, (GLfloat[]) { 0.9f, 0.85f, 0.8f, 0.0f });
	double scale = 0.008;
	double points_width = 0.02;

	double timestep = 0.5;
	double maxspeed = 1;
	int optimal_verlet_steps;

	int div = 1;

	if (seed == 0) {
		// Seed the random
		seed = time(NULL);
	}
	srand(seed);
	if (!data_filled)
		fillData(data, coord, nPoints);

	double kh = KH;
	if (kh == 0)
		kh = compute_kh(nPoints, radius_algorithm) * (xmax - xmin);
	double L = 0.0;
	if (use_verlet == 1) {
		optimal_verlet_steps = compute_optimal_verlet(timestep, maxspeed, kh);
		//printf("Optimal Verlet %d\n", optimal_verlet_steps);
		if (optimal_verlet_steps == -1) { use_verlet = 0; }
		else { L = optimal_verlet_steps * timestep * maxspeed; }
	}

	// create the gridded data structure
	double width = (xmax - xmin) * scale / (kh + L);
	double height = (ymax - ymin) * scale / (kh + L);
	//printf("Width: %f\n", width);
	//printf("Height: %f\n", height);
	//printf("Cells Filled \n");

	// send data to GPU, and receive reference to those data in a points object
	bov_points_t* particles = bov_particles_new(data, nPoints, GL_STATIC_DRAW);

	// setting particles appearance
	bov_points_set_width(particles, 0.02);
	bov_points_set_outline_width(particles, 0.0025);

	// fit particles in a [-0.8 0.8]x[-0.8 0.8] square
	bov_points_scale(particles, (GLfloat[2]) { scale, scale });
	bov_points_set_pos(particles, (GLfloat[2]) { 0.0, -0.1 });

	// we got 0.2 at the top to write something. The screen goes from -1 to 1
	bov_text_t* msg = bov_text_new(
		(GLubyte[]) {
		"Rendering " xstr(nPoints) " particles"
	},
		GL_STATIC_DRAW);
	bov_text_set_pos(msg, (GLfloat[2]) { -0.95, 0.82 });
	bov_text_set_fontsize(msg, 0.1);

	bov_points_t* gridPoint = NULL;
	bov_points_t* pointset = bov_points_new((float[4][2]) {
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
	cell_for_mac* cellArray = NULL;
	double size = (xmax - xmin) / (kh + L);
	double cellLength = (xmax - xmin) / (size);
	use_cells = use_cells && (kh + L) < (xmax - xmin) / 3.0;
	if (use_cells) {
		gridPoint = find_grid_points((kh + L), 100);
		bov_points_set_param(gridPoint, lineParams);
		bov_points_scale(gridPoint, (GLfloat[2]) { scale, scale });
		bov_points_set_pos(gridPoint, (GLfloat[2]) { 0.0, -0.1 });
	}
	bov_points_set_param(pointset, lineParams);
	bov_points_scale(pointset, (GLfloat[2]) { scale, scale });
	bov_points_set_pos(pointset, (GLfloat[2]) { 0.0, -0.1 });

	int cellCounter = 0;
	neighborhood_for_mac** list_of_neighborhood = malloc(iter * sizeof(neighborhood_for_mac*));
	CHECK_MALLOC(list_of_neighborhood);
	int iterations = 0;
	ladfm* loop = malloc(sizeof(ladfm));
	CHECK_MALLOC(loop);
	while (iterations < iter) {
		if (search_neighborhood) {
			neighborhood_for_mac* nh;
			if (use_verlet && (iterations % optimal_verlet_steps))
				nh = neighborhood_new_for_mac(nPoints, list_of_neighborhood[iterations - 1]);
			else
				nh = neighborhood_new_for_mac(nPoints, NULL);

			if (use_cells) {
				cellCounter = 0;
				cellArray = cell_new_for_mac(ceil(size) * ceil(size));
				for (int i = 0; i < nPoints; i++) {
					int cellNumber = ((int)((coord[i][1] - xmin) / (xmax - xmin) * size) * ceil(size) + (int)((coord[i][0] - xmin) / (xmax - xmin) * size));
					if (data[i][1] == xmax)
						cellNumber -= size;
					if (data[i][0] == xmax)
						cellNumber -= 1;
					node_new_for_mac(cellArray, cellNumber, i);
				}
			}
				loop->cells = cellArray;
				loop->use_cells = use_cells;
				loop->use_verlet = use_verlet;
				loop->use_improved_method = use_improved_method;
				loop->div = div;
				loop->size = size;
				loop->nh = nh;
				loop->kh = kh;
				loop->L = L;
				loop->nPoints = nPoints;
				loop->coord = coord;
				loop->data = data;
				loop->scale = scale;
				loop->cellCounter = cellCounter;
				if (use_verlet)
					loop->iterations = iterations % optimal_verlet_steps;
				else
					loop->iterations = 0;
				loop_with_drawing_for_mac(loop);
			while (!bov_window_should_close(window)) {
				bov_particles_draw(window, particles, 0, BOV_TILL_END);
				bov_text_draw(window, msg);
				bov_line_loop_draw(window, pointset, 0, BOV_TILL_END);
				if (use_cells)
					bov_lines_draw(window, gridPoint, 0, BOV_TILL_END);
				bov_window_update(window);
			}
			// Moves the particles around
			if (use_cells)
				cell_delete_for_mac(cellArray, ceil(size) * ceil(size));

			list_of_neighborhood[iterations] = nh;
		}
		bouncyrandomupdate_for_mac(data, coord, 1, -100, 100, -100, 100, 2);
		particles = bov_particles_update(particles, data, NPTS);
		bov_particles_draw(window, particles, 0, BOV_TILL_END);
		bov_text_draw(window, msg);
		bov_line_loop_draw(window, pointset, 0, BOV_TILL_END);
		if (use_cells)
			bov_lines_draw(window, gridPoint, 0, BOV_TILL_END);
		bov_window_update(window);
		iterations++;
	}
	free(loop);

	if (use_cells) {
		bov_points_delete(gridPoint);
	}


	bov_points_delete(pointset);
	bov_text_delete(msg);
	bov_points_delete(particles);
	bov_window_delete(window);

	return list_of_neighborhood;
}

neighborhood_for_mac** create_neighborhood_for_mac(GLfloat(*data)[8], GLfloat(*coord)[2], int iter, int data_filled, int nPoints, int use_cells, int use_improved_method, int radius_algorithm, int use_verlet, unsigned int seed)
{
	//Defining the domain
	int xmax = 100;
	int xmin = -100;
	int ymax = 100;
	int ymin = -100;

	double points_width = 0.02;

	double timestep = 0.5;
	double maxspeed = 1;
	int optimal_verlet_steps;

	if (seed == 0) {
		// Seed the random
		seed = time(NULL);
	}
	srand(seed);
	if (!data_filled)
		fillData_for_mac(data, coord, nPoints);

	double kh = KH;
	if (kh == 0)
		kh = compute_kh_for_mac(nPoints, radius_algorithm) * (xmax - xmin);
	double L = 0.0;
	if (use_verlet == 1) {
		optimal_verlet_steps = compute_optimal_verlet_for_mac(timestep, maxspeed, kh);
		if (optimal_verlet_steps == -1) { use_verlet = 0; }
		else { L = optimal_verlet_steps * timestep * maxspeed; }
	}

	// create the gridded data structure
	double width = (xmax - xmin) / (kh + L);
	double height = (ymax - ymin) / (kh + L);

	cell_for_mac* cellArray = NULL;
	double size = (xmax - xmin) / (kh + L);
	double cellLength = (xmax - xmin) / (size);
	use_cells = use_cells && (kh + L) < (xmax - xmin) / 3.0;
	int cellCounter = 0;
	neighborhood_for_mac** list_of_neighborhood = malloc(iter * sizeof(neighborhood_for_mac*));
	CHECK_MALLOC(list_of_neighborhood);
	int iterations = 0;

	lafm* loop = malloc(sizeof(lafm));
	CHECK_MALLOC(loop);
	while (iterations < iter) {
		neighborhood_for_mac* nh;
		if (use_verlet && (iterations % optimal_verlet_steps))
			nh = neighborhood_new_for_mac(nPoints, list_of_neighborhood[iterations - 1]);
		else
			nh = neighborhood_new_for_mac(nPoints, NULL);

		if (use_cells) {
			cellCounter=0;
			cellArray = cell_new_for_mac(ceil(size) * ceil(size));
			for (int i = 0; i < nPoints; i++) {
				int cellNumber = ((int)((coord[i][1] - xmin) / (xmax - xmin) * size) * ceil(size) + (int)((coord[i][0] - xmin) / (xmax - xmin) * size));
				if (data[i][1] == xmax)
					cellNumber -= size;
				if (data[i][0] == xmax)
					cellNumber -= 1;
				node_new_for_mac(cellArray, cellNumber, i);
			}
		}
			loop->cells = cellArray;
			loop->use_cells = use_cells;
			loop->use_verlet = use_verlet;
			loop->use_improved_method = use_improved_method;
			loop->size = size;
			loop->nh = nh;
			loop->kh = kh;
			loop->L = L;
			loop->nPoints = nPoints;
			loop->coord = coord;
			loop->data = data;
			loop->cellCounter = cellCounter;
			if (use_verlet)
				loop->iterations = iterations % optimal_verlet_steps;
			else
				loop->iterations = 0;
		loop_without_drawing_for_mac(loop);
		// Moves the particles around
		if (use_cells)
			cell_delete_for_mac(cellArray, ceil(size) * ceil(size));
		list_of_neighborhood[iterations] = nh;
		bouncyrandomupdate_for_mac(data, coord, 1, -100, 100, -100, 100, 2);
		iterations++;
	}

	free(loop);
	return list_of_neighborhood;
}

int compare_neighborhood_lists_for_mac(neighborhood_for_mac** nhList_1, neighborhood_for_mac** nhList_2, int iter, int nPoints) {
	int is_equal = 0;
	for (int i = 0; (i < iter) && !is_equal; i++) {
		is_equal = compare_neighborhoods_for_mac(nhList_1[i], nhList_2[i], nPoints);
	}
	return is_equal;
}

int compare_neighborhoods_for_mac(neighborhood_for_mac* nh_1, neighborhood_for_mac* nh_2, int nPoints) {
	for (int i = 0; i < nPoints; i++) {
		if (nh_1[i].nNeighbours != nh_2[i].nNeighbours)
			return 0;
		neighbours_for_mac* current = NULL;
		for (current = nh_1[i].list; current; current = current->next) {
			int is_equal = 0;
			neighbours_for_mac* temp = NULL;
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

void show_equivalence_of_major_algorithm_combination_for_mac(int drawing) {
	GLfloat(*data)[8] = malloc(sizeof(data[0]) * NPTS);
	CHECK_MALLOC(data);
	GLfloat(*coord)[2] = malloc(sizeof(coord[0]) * NPTS);
	CHECK_MALLOC(coord);

	// Seed the random
	time_t seed = time(NULL);
	printf(" %u ", seed);
	srand(seed);

	neighborhood_for_mac** nhList0, ** nhList1, ** nhList2, ** nhList3, ** nhList4, ** nhList5, ** nhList6, ** nhList7, ** nhList8, ** nhList9, ** nhList10, ** nhList11, ** nhList12, ** nhList13, ** nhList14, ** nhList15;
	
	/* Takes way more time and memory, but is supposed to give the same result as whithout drawing*/
	if (drawing) {
		nhList0 = create_neighborhood_with_drawing_for_mac(data, coord, MAX_ITER, 0, NPTS, 1, 1, RADIUS_ALGORITHM, SEARCH_NEIGHBORHOOD, 1, seed);
		nhList1 = create_neighborhood_with_drawing_for_mac(data, coord, MAX_ITER, 0, NPTS, 1, 1, RADIUS_ALGORITHM, SEARCH_NEIGHBORHOOD, 1, seed);
		nhList2 = create_neighborhood_with_drawing_for_mac(data, coord, MAX_ITER, 0, NPTS, 1, 1, RADIUS_ALGORITHM, SEARCH_NEIGHBORHOOD, 0, seed);
		nhList3 = create_neighborhood_with_drawing_for_mac(data, coord, MAX_ITER, 0, NPTS, 1, 1, RADIUS_ALGORITHM, SEARCH_NEIGHBORHOOD, 0, 0, 1, seed);
		nhList4 = create_neighborhood_with_drawing_for_mac(data, coord, MAX_ITER, 0, NPTS, 1, 0, RADIUS_ALGORITHM, SEARCH_NEIGHBORHOOD, 1, seed);
		nhList5 = create_neighborhood_with_drawing_for_mac(data, coord, MAX_ITER, 0, NPTS, 1, 0, RADIUS_ALGORITHM, SEARCH_NEIGHBORHOOD, 1, seed);
		nhList6 = create_neighborhood_with_drawing_for_mac(data, coord, MAX_ITER, 0, NPTS, 1, 0, RADIUS_ALGORITHM, SEARCH_NEIGHBORHOOD, 0, seed);
		nhList7 = create_neighborhood_with_drawing_for_mac(data, coord, MAX_ITER, 0, NPTS, 1, 0, RADIUS_ALGORITHM, SEARCH_NEIGHBORHOOD, 0, seed);
		nhList8 = create_neighborhood_with_drawing_for_mac(data, coord, MAX_ITER, 0, NPTS, 0, 1, RADIUS_ALGORITHM, SEARCH_NEIGHBORHOOD, 1, seed);
		nhList9 = create_neighborhood_with_drawing_for_mac(data, coord, MAX_ITER, 0, NPTS, 0, 1, RADIUS_ALGORITHM, SEARCH_NEIGHBORHOOD, 1, seed);
		nhList10 = create_neighborhood_with_drawing_for_mac(data, coord, MAX_ITER, 0, NPTS, 0, 1, RADIUS_ALGORITHM, SEARCH_NEIGHBORHOOD, 0, seed);
		nhList11 = create_neighborhood_with_drawing_for_mac(data, coord, MAX_ITER, 0, NPTS, 0, 1, RADIUS_ALGORITHM, SEARCH_NEIGHBORHOOD, 0, seed);
		nhList12 = create_neighborhood_with_drawing_for_mac(data, coord, MAX_ITER, 0, NPTS, 0, 0, RADIUS_ALGORITHM, SEARCH_NEIGHBORHOOD, 1, seed);
		nhList13 = create_neighborhood_with_drawing_for_mac(data, coord, MAX_ITER, 0, NPTS, 0, 0, RADIUS_ALGORITHM, SEARCH_NEIGHBORHOOD, 1, seed);
		nhList14 = create_neighborhood_with_drawing_for_mac(data, coord, MAX_ITER, 0, NPTS, 0, 0, RADIUS_ALGORITHM, SEARCH_NEIGHBORHOOD, 0, seed);
		nhList15 = create_neighborhood_with_drawing_for_mac(data, coord, MAX_ITER, 0, NPTS, 0, 0, RADIUS_ALGORITHM, SEARCH_NEIGHBORHOOD, 0, seed);
	}
	else {
		nhList0 = create_neighborhood_for_mac(data, coord, MAX_ITER, 0, NPTS, 1, 1, RADIUS_ALGORITHM, 1, seed);
		nhList1 = create_neighborhood_for_mac(data, coord, MAX_ITER, 0, NPTS, 1, 1, RADIUS_ALGORITHM, 1, seed);
		nhList2 = create_neighborhood_for_mac(data, coord, MAX_ITER, 0, NPTS, 1, 1, RADIUS_ALGORITHM, 0, seed);
		nhList3 = create_neighborhood_for_mac(data, coord, MAX_ITER, 0, NPTS, 1, 1, RADIUS_ALGORITHM, 0, seed);
		nhList4 = create_neighborhood_for_mac(data, coord, MAX_ITER, 0, NPTS, 1, 0, RADIUS_ALGORITHM, 1, seed);
		nhList5 = create_neighborhood_for_mac(data, coord, MAX_ITER, 0, NPTS, 1, 0, RADIUS_ALGORITHM, 1, seed);
		nhList6 = create_neighborhood_for_mac(data, coord, MAX_ITER, 0, NPTS, 1, 0, RADIUS_ALGORITHM, 0, seed);
		nhList7 = create_neighborhood_for_mac(data, coord, MAX_ITER, 0, NPTS, 1, 0, RADIUS_ALGORITHM, 0, seed);
		nhList8 = create_neighborhood_for_mac(data, coord, MAX_ITER, 0, NPTS, 0, 1, RADIUS_ALGORITHM, 1, seed);
		nhList9 = create_neighborhood_for_mac(data, coord, MAX_ITER, 0, NPTS, 0, 1, RADIUS_ALGORITHM, 1, 0, 1, seed);
		nhList10 = create_neighborhood_for_mac(data, coord, MAX_ITER, 0, NPTS, 0, 1, RADIUS_ALGORITHM, 0, seed);
		nhList11 = create_neighborhood_for_mac(data, coord, MAX_ITER, 0, NPTS, 0, 1, RADIUS_ALGORITHM, 0, seed);
		nhList12 = create_neighborhood_for_mac(data, coord, MAX_ITER, 0, NPTS, 0, 0, RADIUS_ALGORITHM, 1, seed);
		nhList13 = create_neighborhood_for_mac(data, coord, MAX_ITER, 0, NPTS, 0, 0, RADIUS_ALGORITHM, 1, seed);
		nhList14 = create_neighborhood_for_mac(data, coord, MAX_ITER, 0, NPTS, 0, 0, RADIUS_ALGORITHM, 0, seed);
		nhList15 = create_neighborhood_for_mac(data, coord, MAX_ITER, 0, NPTS, 0, 0, RADIUS_ALGORITHM, 0, seed);
	}

	printf("\n0==1 ?");
	printf("\n  %i  \n", compare_neighborhood_lists_for_mac(nhList0, nhList1, MAX_ITER, NPTS));
	printf("\n1==2 ?");
	printf("\n  %i  \n", compare_neighborhood_lists_for_mac(nhList1, nhList2, MAX_ITER, NPTS));
	printf("\n2==3 ?");
	printf("\n  %i  \n", compare_neighborhood_lists_for_mac(nhList2, nhList3, MAX_ITER, NPTS));
	printf("\n3==4 ?");
	printf("\n  %i  \n", compare_neighborhood_lists_for_mac(nhList3, nhList4, MAX_ITER, NPTS));
	printf("\n4==5 ?");
	printf("\n  %i  \n", compare_neighborhood_lists_for_mac(nhList4, nhList5, MAX_ITER, NPTS));
	printf("\n5==6 ?");
	printf("\n  %i  \n", compare_neighborhood_lists_for_mac(nhList5, nhList6, MAX_ITER, NPTS));
	printf("\n6==7 ?");
	printf("\n  %i  \n", compare_neighborhood_lists_for_mac(nhList6, nhList7, MAX_ITER, NPTS));
	printf("\n7==8 ?");
	printf("\n  %i  \n", compare_neighborhood_lists_for_mac(nhList7, nhList8, MAX_ITER, NPTS));
	printf("\n8==9 ?");
	printf("\n  %i  \n", compare_neighborhood_lists_for_mac(nhList8, nhList9, MAX_ITER, NPTS));
	printf("\n9==10 ?");
	printf("\n  %i  \n", compare_neighborhood_lists_for_mac(nhList9, nhList10, MAX_ITER, NPTS));
	printf("\n10==11 ?");
	printf("\n  %i  \n", compare_neighborhood_lists_for_mac(nhList10, nhList11, MAX_ITER, NPTS));
	printf("\n11==12 ?");
	printf("\n  %i  \n", compare_neighborhood_lists_for_mac(nhList11, nhList12, MAX_ITER, NPTS));
	printf("\n12==13 ?");
	printf("\n  %i  \n", compare_neighborhood_lists_for_mac(nhList12, nhList13, MAX_ITER, NPTS));
	printf("\n13==14 ?");
	printf("\n  %i  \n", compare_neighborhood_lists_for_mac(nhList13, nhList14, MAX_ITER, NPTS));
	printf("\n14==15 ?");
	printf("\n  %i  \n", compare_neighborhood_lists_for_mac(nhList14, nhList15, MAX_ITER, NPTS));

	for (int i = 0; i < MAX_ITER; i++) {
		neighborhood_delete_for_mac(nhList0[i], NPTS);
	}
	for (int i = 0; i < MAX_ITER; i++) {
		neighborhood_delete_for_mac(nhList1[i], NPTS);
	}
	for (int i = 0; i < MAX_ITER; i++) {
		neighborhood_delete_for_mac(nhList2[i], NPTS);
	}
	for (int i = 0; i < MAX_ITER; i++) {
		neighborhood_delete_for_mac(nhList3[i], NPTS);
	}
	for (int i = 0; i < MAX_ITER; i++) {
		neighborhood_delete_for_mac(nhList4[i], NPTS);
	}
	for (int i = 0; i < MAX_ITER; i++) {
		neighborhood_delete_for_mac(nhList5[i], NPTS);
	}
	for (int i = 0; i < MAX_ITER; i++) {
		neighborhood_delete_for_mac(nhList6[i], NPTS);
	}
	for (int i = 0; i < MAX_ITER; i++) {
		neighborhood_delete_for_mac(nhList7[i], NPTS);
	}
	for (int i = 0; i < MAX_ITER; i++) {
		neighborhood_delete_for_mac(nhList8[i], NPTS);
	}
	for (int i = 0; i < MAX_ITER; i++) {
		neighborhood_delete_for_mac(nhList9[i], NPTS);
	}
	for (int i = 0; i < MAX_ITER; i++) {
		neighborhood_delete_for_mac(nhList10[i], NPTS);
	}
	for (int i = 0; i < MAX_ITER; i++) {
		neighborhood_delete_for_mac(nhList11[i], NPTS);
	}
	for (int i = 0; i < MAX_ITER; i++) {
		neighborhood_delete_for_mac(nhList12[i], NPTS);
	}
	for (int i = 0; i < MAX_ITER; i++) {
		neighborhood_delete_for_mac(nhList13[i], NPTS);
	}
	for (int i = 0; i < MAX_ITER; i++) {
		neighborhood_delete_for_mac(nhList14[i], NPTS);
	}
	for (int i = 0; i < MAX_ITER; i++) {
		neighborhood_delete_for_mac(nhList15[i], NPTS);
	}
	free(nhList0);
	free(nhList1);
	free(nhList2);
	free(nhList3);
	free(nhList4);
	free(nhList5);
	free(nhList6);
	free(nhList7);
	free(nhList8);
	free(nhList9);
	free(nhList10);
	free(nhList11);
	free(nhList12);
	free(nhList13);
	free(nhList14);
	free(nhList15);
	free(coord);
	free(data);
}