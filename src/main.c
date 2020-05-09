#include "neighborhood_search.h"
#include "kernel.h"
#include "conservation_equations.h"
#include <math.h> 
#include <time.h> 
#include <stdio.h> 
/*
int NPTS = 1600;//225;
double MASS = 5;
double VELOCITY = 1;
double INITDENSITY = 0.035;// 0.05;
double ID = 0.0349;
int IS_CIRCULAR = 0;
int TEMPERATURE = 0;
int USE_EULER = 0;
*/
int NPTS = 225;
double MASS = 5;
double VELOCITY = 1000;
double INITDENSITY = 0.5;
double ID = 0.005;
int IS_CIRCULAR = 1;
int TEMPERATURE = 0;
int USE_EULER = 0;

// for v, a floating point value between 0 and 1, this function fills color with
// the improved jet colormap color corresponding to v
static void circular_velocity_to_colormap(float v, float color[3])
{
	color[0] = 1.5 - 4.0 * fabs(v / 0.2 / VELOCITY - 0.25);
	color[1] = 1.5 - 4.0 * fabs(v / 0.2 / VELOCITY + 0.0);
	color[2] = 1.5 - 4.0 * fabs(v / 0.2 / VELOCITY + 0.25);
}
static void square_velocity_to_colormap(float v, float color[3])
{
	color[0] = 1.5 - 4.0 * fabs(v / 0.8 / VELOCITY - 0.75);
	color[1] = 1.5 - 4.0 * fabs(v / 0.8 / VELOCITY - 0.5);
	color[2] = 1.5 - 4.0 * fabs(v / 0.8 / VELOCITY - 0.25);
}
static void temperature_to_colormap(float v, float color[3])
{
	color[0] = 1.5 - 4.0 * fabs((v-323.15)/ 100 - 0.25);
	color[1] = 1.5 - 4.0 * fabs((v - 323.15) / 100 - 0.0);
	color[2] = 1.5 - 4.0 * fabs((v - 323.15) / 100 + 0.25);
}

// function to fill the data table of the nPoints particles positions, speeds, colors and transparency and the coord table with the nPoints particles positions used to draw;
// data[i][0] == coord[i][0] && data[i][1] == coord[i][1]
void circular_fillData(GLfloat(*data)[8], double out_r, double in_r)
{
	float rmax = (out_r - 2.5) * sqrtf(2.0f);
	for (int i = 0; i < NPTS; i++) {
		double theta = rand() * 2 * M_PI / RAND_MAX;
		double radius = rand() * (out_r - in_r -5) / RAND_MAX + in_r+2.5;
		data[i][0] = radius*cos(theta); // x (rand between 2 circles)
		data[i][1] = radius * sin(theta); // y (rand between 2 circles)
		double r = sqrt(data[i][0] * data[i][0] + data[i][1] * data[i][1]);
		data[i][2] = 0;// rand() * 2.0 / RAND_MAX - 1.0;; //Random starting speed
		data[i][3] = 0; // rand() * 2.0 / RAND_MAX - 1.0; //Random starting speed
		circular_velocity_to_colormap(0, &data[i][4]); // fill color
		data[i][7] = 0.8f; // transparency
	}
}

void random_square_fillData(GLfloat(*data)[8], double half_length)
{
	for (int i = 0; i < NPTS; i++) {
		data[i][0] = rand() * 2 * half_length / RAND_MAX - half_length; // x (rand between -100 and 100)
		data[i][1] = rand() * 2 * half_length / RAND_MAX - half_length; // y (rand between -100 and 100)
		data[i][2] = 0;
		data[i][3] = 0;
		square_velocity_to_colormap(0, &data[i][4]); // fill color
		data[i][7] = 0.8f; // transparency
	}
}

void square_fillData(GLfloat(*data)[8], double half_length)
{
	int height = ceil(sqrt(NPTS));
	int width = ceil(sqrt(NPTS));
	for (int i = 0; i < height; i++) {
		for (int j = 0; j < width; j++) {
			data[i * height + j][0] = -half_length + j * 2 * half_length / (width - 1); // x (rand between -100 and 100)
			data[i * height + j][1] = half_length - i * 2 * half_length / (height - 1); // y (rand between -100 and 100)
			data[i * height + j][2] = 0;
			data[i * height + j][3] = 0;
			square_velocity_to_colormap(0, &data[i * height + j][4]); // fill color
			data[i * height + j][7] = 0.8f; // transparency
		}
	}
}

void print_points(GLfloat(*data)[4]) {
	for (int i = 0; i < NPTS; i++)
		printf("%i %f %f\n", i + 1, data[i][2], data[i][3]);
}

void  circular_boundary_particles_create(int outside_stages, int inside_stages, GLfloat(*boundary)[8], double out_r, double in_r, double a_c) {
	double r = out_r+2.5;
	int current = 0;
	for (int i = 0; i < outside_stages; i++) {
		int j_max = (int)(2.0 * M_PI * (r + (double)i * 5.0) / 5.0);
		double v = sqrt(a_c * r);
		for (int j = 0; j < j_max; j++) {
			boundary[current + j][0] = (r + (double)i * 5.0) * cos(2 * M_PI * j / j_max);
			boundary[current + j][1] = (r + (double)i * 5.0) * sin(2 * M_PI * j / j_max);
			boundary[current + j][2] = -sin(2 * M_PI * j / j_max)*v;
			boundary[current + j][3] = cos(2 * M_PI * j / j_max)*v;
			boundary[current + j][4] = 1;
			boundary[current + j][5] = 1;
			boundary[current + j][6] = 1;
			boundary[current + j][7] = 1.0;
		}
		current += j_max;
	}
	r = in_r - 2.5;
	for (int i = 0; i < inside_stages; i++) {
		int j_max = (int)(2.0 * M_PI * (r - (double)i * 5.0) / 5.0);
		double v = -sqrt(50*a_c * r);
		for (int j = 0; j < j_max; j++) {
			boundary[current + j][0] = (r - (double)i * 5.0) * cos(2 * M_PI * j / j_max);
			boundary[current + j][1] = (r - (double)i * 5.0) * sin(2 * M_PI * j / j_max);
			boundary[current + j][2] = -sin(2 * M_PI * j / j_max)*v;
			boundary[current + j][3] = cos(2 * M_PI * j / j_max)*v;
			boundary[current + j][4] = 1;
			boundary[current + j][5] = 1;
			boundary[current + j][6] = 1;
			boundary[current + j][7] = 1.0;
		}
		current += j_max;
	}
}

void  square_boundary_particles_create(int n_stages, GLfloat(*boundary)[8], double half_length) {
	int current = 0;
	for (int i = 0; i < n_stages; i++) {
		int j_max = (int)(4 * (n_stages + half_length / 5));
		for (int j = 0; j < j_max; j = j + 2) {
			boundary[current + j][0] = 2.5 + 5 * j / 2 - (n_stages * 5 + half_length);
			boundary[current + j][1] = half_length + 2.5 + i * 5;
			boundary[current + j][2] = VELOCITY;
			boundary[current + j][3] = 0;
			boundary[current + j][4] = 1;
			boundary[current + j][5] = 1;
			boundary[current + j][6] = 1;
			boundary[current + j][7] = 1.0;
			boundary[current + j + 1][0] = 2.5 + 5 * j / 2 - (n_stages * 5 + half_length);
			boundary[current + j + 1][1] = -(half_length + 2.5 + i * 5);
			boundary[current + j + 1][2] = 0;//-VELOCITY;
			boundary[current + j + 1][3] = 0;
			boundary[current + j + 1][4] = 1;
			boundary[current + j + 1][5] = 1;
			boundary[current + j + 1][6] = 1;
			boundary[current + j + 1][7] = 1.0;
		}
		current += j_max;
	}
	for (int i = 0; i < n_stages; i++) {
		int j_max = (int)(4 * (half_length / 5));
		for (int j = 0; j < j_max; j = j + 2) {
			boundary[current + j][0] = half_length + 2.5 + i * 5;
			boundary[current + j][1] = 2.5 + 5 * j / 2 - half_length;
			boundary[current + j][2] = 0;
			boundary[current + j][3] = 0;// -VELOCITY;
			boundary[current + j][4] = 1;
			boundary[current + j][5] = 1;
			boundary[current + j][6] = 1;
			boundary[current + j][7] = 1.0;
			boundary[current + j + 1][0] = -(half_length + 2.5 + i * 5);
			boundary[current + j + 1][1] = 2.5 + 5 * j / 2 - half_length;
			boundary[current + j + 1][2] = 0;
			boundary[current + j + 1][3] = 0;// VELOCITY;
			boundary[current + j + 1][4] = 1;
			boundary[current + j + 1][5] = 1;
			boundary[current + j + 1][6] = 1;
			boundary[current + j + 1][7] = 1.0;
		}
		current += j_max;
	}
}

void  temperature_boundary_particles_create(int n_stages, GLfloat(*boundary)[8], double half_length) {
	int current = 0;
	for (int i = 0; i < n_stages; i++) {
		int j_max = (int)(4 * (n_stages + half_length / 5));
		for (int j = 0; j < j_max; j = j + 2) {
			boundary[current + j][0] = 2.5 + 5 * j / 2 - (n_stages * 5 + half_length);
			boundary[current + j][1] = half_length + 2.5 + i * 5;
			boundary[current + j][2] = 0;
			boundary[current + j][3] = 0;
			boundary[current + j][4] = 1;
			boundary[current + j][5] = 1;
			boundary[current + j][6] = 1;
			boundary[current + j][7] = 1.0;
			boundary[current + j + 1][0] = 2.5 + 5 * j / 2 - (n_stages * 5 + half_length);
			boundary[current + j + 1][1] = -(half_length + 2.5 + i * 5);
			boundary[current + j + 1][2] = 0;//-VELOCITY;
			boundary[current + j + 1][3] = 0;
			boundary[current + j + 1][4] = 1;
			boundary[current + j + 1][5] = 1;
			boundary[current + j + 1][6] = 1;
			boundary[current + j + 1][7] = 1.0;
		}
		current += j_max;
	}
	for (int i = 0; i < n_stages; i++) {
		int j_max = (int)(4 * (half_length / 5));
		for (int j = 0; j < j_max; j = j + 2) {
			boundary[current + j][0] = half_length + 2.5 + i * 5;
			boundary[current + j][1] = 2.5 + 5 * j / 2 - half_length;
			boundary[current + j][2] = 0;
			boundary[current + j][3] = 0;// -VELOCITY;
			boundary[current + j][4] = 1;
			boundary[current + j][5] = 1;
			boundary[current + j][6] = 1;
			boundary[current + j][7] = 1.0;
			boundary[current + j + 1][0] = -(half_length + 2.5 + i * 5);
			boundary[current + j + 1][1] = 2.5 + 5 * j / 2 - half_length;
			boundary[current + j + 1][2] = 0;
			boundary[current + j + 1][3] = 0;// VELOCITY;
			boundary[current + j + 1][4] = 1;
			boundary[current + j + 1][5] = 1;
			boundary[current + j + 1][6] = 1;
			boundary[current + j + 1][7] = 1.0;
		}
		current += j_max;
	}
}

void initial_condition(GLfloat(*sup_data)[3], double volume, int total) {
	printf("%f\n", MASS / volume);
	printf("%f\n", MASS / NPTS);
	for (int i = 0; i < NPTS; i++) {
		sup_data[i][0] = INITDENSITY;// MASS / 1 / volume;
		sup_data[i][1] = 0.05 * volume / NPTS;
		sup_data[i][2] = 285;
	}
	for (int i = NPTS; i < (NPTS + total); i++) {
		sup_data[i][0] = INITDENSITY;// MASS / 1 / volume;
		sup_data[i][1] = 0.05 * volume / NPTS;
		sup_data[i][2] = 335;
	}
}

void temperature_initial_condition(GLfloat(*data)[8],GLfloat(*sup_data)[3], double volume, int total, int half_length) {
	printf("%f\n", MASS / volume);
	printf("%f\n", MASS / NPTS);
	for (int i = 0; i < NPTS; i++) {
		sup_data[i][0] = INITDENSITY;// MASS / 1 / volume;
		sup_data[i][1] = 0.05 * volume / NPTS;
		sup_data[i][2] = 273.15;
	}
	for (int i = NPTS; i < (NPTS + total); i++) {
		sup_data[i][0] = INITDENSITY;// MASS / 1 / volume;
		sup_data[i][1] = 0.05 * volume / NPTS;
		if(data[i][1]<-half_length)
			sup_data[i][2] = 373.15;
		else
			sup_data[i][2] = 273.15;
	}
}

void singular_point_lid(GLfloat(*data)[8], GLfloat(*sup_data)[3], int total, int half_length) {
	for (int i = NPTS; i < total + NPTS; i++) {
		if (data[i][0] < (-half_length) && data[i][1] < (-half_length))
			sup_data[i][0] = ID;
		/*else if (data[i][0] < (-half_length) && data[i][1] > (half_length))
			sup_data[i][0] = ID;
		else if (data[i][0] > (half_length) && data[i][1] < (-half_length))
			sup_data[i][0] = ID;
		else if (data[i][0] > (half_length) && data[i][1] > (half_length))
			sup_data[i][0] = ID;*/
	}
}

int main()
{
	bov_window_t* window = bov_window_new(1024, 780, "ANM Project: SPH");
	bov_window_set_color(window, (GLfloat[]) { 0.9f, 0.85f, 0.8f, 0.0f });
	//print_points(data);
	double timestep = 0.05;
	double maxspeed = 1;
	neighborhood_options* options = neighborhood_options_init(timestep, maxspeed);
	neighborhood* nh = options->nh;
	int number_of_iterations = 10000;
	double radius_min, radius_max, a_c, half_length;
	int total, out, outside_stages,inside_stages,n_stages;
	if (IS_CIRCULAR) {
		radius_min = 20;
		radius_max = 100.0 - 20;

		double r_max = radius_max + 2.5;
		a_c = VELOCITY / r_max;
		total = 0;
		outside_stages = floor(options->kh / 5);
		for (int i = 0; i < outside_stages; i++)
			total += (int)(2.0 * M_PI * (r_max + (double)i * 5.0) / 5.0);
		out = total;
		double r_min = radius_min - 2.5;
		inside_stages = floor(fmin(options->kh, radius_min) / 5);
		for (int i = 0; i < inside_stages; i++)
			total += (int)(2.0 * M_PI * (r_min - (double)i * 5.0) / 5.0);
	}
	else {
		half_length = 100;
		//n_stages = options->kh / 5;
		n_stages = 30 / 5;
		total = 4 * n_stages * (n_stages + (int)(2 * half_length / 5));
	}
	GLfloat(*data)[8] = malloc(sizeof(data[0]) * (NPTS + total));
	CHECK_MALLOC(data);
	GLfloat(*sup_data)[3] = malloc(sizeof(sup_data[0]) * (NPTS+total));
	CHECK_MALLOC(sup_data);
	GLfloat(*new_data)[6] = calloc(NPTS, sizeof(new_data[0]));
	CHECK_MALLOC(new_data);
	for (int i = 0; i < NPTS; i++)
		for (int j = 0; j < 4; j++)
			new_data[i][j] = 0.0;
	//print_points(new_data);
	// Seed the random
	time_t seed = time(NULL);
	//printf(" %u \n", seed);
	srand(seed);

	if (IS_CIRCULAR) {
		circular_fillData(data, radius_max, radius_min);
		circular_boundary_particles_create(outside_stages,inside_stages, &data[NPTS], radius_max,radius_min,a_c);
		initial_condition(sup_data,M_PI*(pow(radius_max,2)- pow(radius_min,2)),total);
	}
	else {
		square_fillData(data, half_length-5);
		if (TEMPERATURE) {
			temperature_boundary_particles_create(n_stages, &data[NPTS], half_length);
			temperature_initial_condition(data,sup_data, pow(half_length, 2), total,half_length);
		}
		else {
			square_boundary_particles_create(n_stages, &data[NPTS], half_length);
			initial_condition(sup_data, pow(half_length, 2), total);
			singular_point_lid(data, sup_data, total, half_length);
		}
	}
	GLfloat(*inner_boundary)[2];
	GLfloat(*outer_boundary)[2];
	GLfloat(*boundary)[2];
	if (IS_CIRCULAR) {
		inner_boundary = malloc(sizeof(inner_boundary[0]) * 100);
		CHECK_MALLOC(inner_boundary);
		for (int i = 0; i < 100; i++) {
			inner_boundary[i][0] = radius_min * cos(i * 2 * M_PI / 100);
			inner_boundary[i][1] = radius_min * sin(i * 2 * M_PI / 100);
		}
		outer_boundary = malloc(sizeof(outer_boundary[0]) * 300);
		CHECK_MALLOC(outer_boundary);
		for (int i = 0; i < 300; i++) {
			outer_boundary[i][0] = radius_max * cos(i * 2 * M_PI / 300);
			outer_boundary[i][1] = radius_max * sin(i * 2 * M_PI / 300);
		}
	}
	else {
		boundary = malloc(sizeof(boundary[0]) * 4);
		CHECK_MALLOC(boundary);
		boundary[0][0] = half_length;
		boundary[0][1] = half_length;
		boundary[1][0] = half_length;
		boundary[1][1] = -half_length;
		boundary[2][0] = -half_length;
		boundary[2][1] = -half_length;
		boundary[3][0] = -half_length;
		boundary[3][1] = half_length;
	}

	/* send data to GPU, and receive reference to those data in a points object */
	bov_points_t* particles = bov_particles_new(data, NPTS+total, GL_STATIC_DRAW);
	bov_points_t* inner_particles, * outer_particles, * boundary_particles;
	if (IS_CIRCULAR) {
		inner_particles = bov_points_new(inner_boundary, 100, GL_STATIC_DRAW);
		outer_particles = bov_points_new(outer_boundary, 300, GL_STATIC_DRAW);
	}
	else
		boundary_particles = bov_points_new(boundary, 4, GL_STATIC_DRAW);

	/* setting particles appearance */
	bov_points_set_width(particles, 0.02);
	bov_points_set_outline_width(particles, 0.0025);
	if (IS_CIRCULAR) {
		bov_points_set_width(inner_particles, 0.005);
		bov_points_set_outline_width(inner_particles, 0.0025);
		bov_points_set_width(outer_particles, 0.005);
		bov_points_set_outline_width(outer_particles, 0.0025);
	}
	else {
		bov_points_set_width(boundary_particles, 0.005);
		bov_points_set_outline_width(boundary_particles, 0.0025);
	}
	
	/* fit particles in a [-0.8 0.8]x[-0.8 0.8] square */
	bov_points_scale(particles, (GLfloat[2]) { 0.008, 0.008 });
	bov_points_set_pos(particles, (GLfloat[2]) { 0.0, -0.1 });
	if (IS_CIRCULAR) {
		bov_points_scale(inner_particles, (GLfloat[2]) { 0.008, 0.008 });
		bov_points_set_pos(inner_particles, (GLfloat[2]) { 0.0, -0.1 });
		bov_points_scale(outer_particles, (GLfloat[2]) { 0.008, 0.008 });
		bov_points_set_pos(outer_particles, (GLfloat[2]) { 0.0, -0.1 });
	}
	else {
		bov_points_scale(boundary_particles, (GLfloat[2]) { 0.008, 0.008 });
		bov_points_set_pos(boundary_particles, (GLfloat[2]) { 0.0, -0.1 });
	}
	char npoints[10];
	sprintf(npoints, "%i", NPTS);
	int u = NPTS;
	/* we got 0.2 at the top to write something. The screen goes from -1 to 1 */
	bov_text_t* msg = bov_text_new(
		(GLubyte[]) {
		"Rendering " xstr(200) " particles"
	},
		GL_STATIC_DRAW);
	bov_text_set_pos(msg, (GLfloat[2]) { -0.95, 0.82 });
	bov_text_set_fontsize(msg, 0.1);

	clock_t init_t = clock();
	clock_t prev_t = init_t;
	double final_t = ((double)init_t)+300* CLOCKS_PER_SEC;
	while((double)clock()<=final_t){
		timestep = ((double)(init_t - prev_t))/CLOCKS_PER_SEC;
		if (timestep != 0) {
			if(TEMPERATURE){
				if (USE_EULER)
					update_euler(data, sup_data, new_data, nh, timestep, options->kh, grad_w_cubic, 3);
				else {
					update_predictor(data, sup_data, new_data, nh, timestep, options->kh, grad_w_cubic, 3);
					timestep /= 2;
				}
				for (int i = 0; i < NPTS; i++) {
					sup_data[i][2] = new_data[i][1];
					temperature_to_colormap(sup_data[i][2], &data[i][4]);
				}
			}
			else {
				update_euler(data, sup_data, new_data, nh, timestep, options->kh, w_cubic, 1);
				if (USE_EULER)
					update_euler(data, sup_data, new_data, nh, timestep, options->kh, grad_w_cubic, 2);
				else {
					update_predictor(data, sup_data, new_data, nh, timestep, options->kh, grad_w_cubic, 2);
					timestep /= 2;
				}
				for (int i = 0; i < NPTS; i++) {
					//printf("%f %f %f %f\n", data[i][0], data[i][1], new_data[i][4], new_data[i][5]);
					sup_data[i][0] = new_data[i][0];
					//printf("%i %f\n", i, new_data[i][0]);
					data[i][2] = new_data[i][2] + new_data[i][4] / sup_data[i][1] * timestep;
					data[i][3] = new_data[i][3] + new_data[i][5] / sup_data[i][1] * timestep;
					//if (sqrt(pow(data[i][2], 2) + pow(data[i][3], 2)) > VELOCITY) {
					//	data[i][2] *= VELOCITY / sqrt(pow(data[i][2], 2) + pow(data[i][3], 2));
					//	data[i][3] *= VELOCITY / sqrt(pow(data[i][2], 2) + pow(data[i][3], 2));
					//}
					data[i][0] += timestep * (new_data[i][2] + new_data[i][4] / sup_data[i][1] * timestep);
					data[i][1] += timestep * (new_data[i][3] + new_data[i][5] / sup_data[i][1] * timestep);
					if (IS_CIRCULAR) {
						double r = sqrt(pow(data[i][0], 2) + pow(data[i][1], 2));
						double cos_theta = data[i][0] / r;
						double sin_theta = data[i][1] / r;
						double tan_theta_v = data[i][3] / data[i][2];
						double v = sqrt(pow(data[i][2], 2) + pow(data[i][3], 2));
						double cos_theta_v = data[i][2] / v;
						double sin_theta_v = data[i][3] / v;
						if (cos_theta < 0 && cos_theta_v > 0)
							circular_velocity_to_colormap(v * sin(atan(sin_theta / cos_theta) - atan(sin_theta_v / cos_theta_v) + M_PI), &data[i][4]);
						else if (cos_theta > 0 && cos_theta_v < 0)
							circular_velocity_to_colormap(v * sin(atan(sin_theta / cos_theta) - atan(sin_theta_v / cos_theta_v) - M_PI), &data[i][4]);
						else
							circular_velocity_to_colormap(v * sin(atan(sin_theta / cos_theta) - atan(sin_theta_v / cos_theta_v)), &data[i][4]);
					}
					else
						square_velocity_to_colormap(sqrt(pow(data[i][2], 2) + pow(data[i][3], 2)), &data[i][4]);
				}
			}
			if (!USE_EULER)
				timestep *= 2;
		}
		neighborhood_update(options, nh, data, 0,total);
		/*if (IS_CIRCULAR) {
			for (int i = 0; i < out; i++) {
				double r = sqrt(pow(data[i + NPTS][0], 2) + pow(data[i + NPTS][1], 2));
				double v = sqrt(a_c * r);
				double cos_theta = data[i + NPTS][0] / r;
				double sin_theta = data[i + NPTS][1] / r;
				data[i + NPTS][0] += -sin_theta * v * timestep + pow(timestep, 2) * cos_theta / (2 * r);
				data[i + NPTS][1] += cos_theta * v * timestep + pow(timestep, 2) * sin_theta / (2 * r);
				double r_new = sqrt(pow(data[i + NPTS][0], 2) + pow(data[i + NPTS][1], 2));
				data[i + NPTS][0] *= r / r_new;
				data[i + NPTS][1] *= r / r_new;
				data[i + NPTS][2] += -sin_theta * v;
				data[i + NPTS][3] += cos_theta * v;
			}
			for (int i = out; i < total; i++) {
				double r = sqrt(pow(data[i + NPTS][0], 2) + pow(data[i + NPTS][1], 2));
				double v = -sqrt(0.3 * a_c * r);
				double cos_theta = data[i + NPTS][0] / r;
				double sin_theta = data[i + NPTS][1] / r;
				data[i + NPTS][0] += -sin_theta * v * timestep + pow(timestep, 2) * cos_theta / (2 * r);
				data[i + NPTS][1] += cos_theta * v * timestep + pow(timestep, 2) * sin_theta / (2 * r);
				double r_new = sqrt(pow(data[i + NPTS][0], 2) + pow(data[i + NPTS][1], 2));
				data[i + NPTS][0] *= r / r_new;
				data[i + NPTS][1] *= r / r_new;
				data[i + NPTS][2] += -sin_theta * v;
				data[i + NPTS][3] += cos_theta * v;
			}
		}
		else {
			for (int i = 0; i < total; i++) {
				data[i + NPTS][0] += data[i + NPTS][2] * timestep;
				if (data[i + NPTS][0] > (half_length + 5 * n_stages))
					data[i + NPTS][0] -= 2 * (half_length + 5 * n_stages);
			}
		}*/
		bov_particles_update(particles, data, NPTS+total);
		bov_particles_draw(window, particles, 0, BOV_TILL_END);
		if (IS_CIRCULAR) {
			bov_line_loop_draw(window, inner_particles, 0, BOV_TILL_END);
			bov_line_loop_draw(window, outer_particles, 0, BOV_TILL_END);
		}
		else
			bov_line_loop_draw(window, boundary_particles, 0, BOV_TILL_END); 

		bov_text_draw(window, msg);

		// In your actual project, don't wait for events => bov_window_update(window)
		bov_window_update(window);
		prev_t = init_t;
		init_t = clock();
	}
	
	while (!bov_window_should_close(window)) {
		bov_particles_draw(window, particles, 0, BOV_TILL_END);
		/*bov_line_loop_draw(window, inner_particles, 0, BOV_TILL_END);
		bov_line_loop_draw(window, outer_particles, 0, BOV_TILL_END);*/
		bov_line_loop_draw(window, boundary_particles, 0, BOV_TILL_END);
		// bov_points_draw(window, particles, 0, BOV_TILL_END);

		bov_text_draw(window, msg);

		// In your actual project, don't wait for events => bov_window_update(window)
		bov_window_update(window);
	}
	neighborhood_options_delete(options,nh);

	bov_text_delete(msg);
	bov_points_delete(particles);
	bov_window_delete(window);

	free(data);
	return EXIT_SUCCESS;
}
