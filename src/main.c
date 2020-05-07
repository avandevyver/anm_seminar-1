#include "neighborhood_search.h"
#include "kernel.h"
#include "conservation_equations.h"
#include <math.h> 
#include <time.h> 
#include <stdio.h> 

int NPTS = 200;
double MASS = 0.5;
double VELOCITY = 0.5;
// for v, a floating point value between 0 and 1, this function fills color with
// the improved jet colormap color corresponding to v
static void velocity_to_colormap(float v, float color[3])
{
	//float v1 = 3.5 * (v - 0.7);
	//float v2 = 1.25 * v;
	//float v3 = fminf(0.5, v) * 2.0;

	//color[0] = -v1 * v1 + 1.0f;
	//color[1] = 6.0f * v2 * v2 * (1.0f - v2);
	//color[2] = 5.5f * v3 * (1.0f - v3) * (1.0f - v3);

	// alternative: classical jet colormap
	 color[0] = 1.5 - 4.0 * fabs(v/2.5/VELOCITY - 0.25);
	 color[1] = 1.5 - 4.0 * fabs(v /2.5/VELOCITY + 0.0 );
	 color[2] = 1.5 - 4.0 * fabs(v /2.5/ VELOCITY + 0.25);
}

// function to fill the data table of the nPoints particles positions, speeds, colors and transparency and the coord table with the nPoints particles positions used to draw;
// data[i][0] == coord[i][0] && data[i][1] == coord[i][1]
void fillData(GLfloat(*data)[8], double out_r, double in_r)
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
		velocity_to_colormap(0, &data[i][4]); // fill color
		data[i][7] = 0.8f; // transparency
	}
}

void print_points(GLfloat(*data)[4]) {
	for (int i = 0; i < NPTS; i++)
		printf("%i %f %f\n", i + 1, data[i][2], data[i][3]);
}

void  boundary_particles_create(int outside_stages, int inside_stages, GLfloat(*boundary)[8], double out_r, double in_r, double a_c) {
	double r = out_r+2.5;
	int current = 0;
	for (int i = 0; i < outside_stages; i++) {
		int j_max = (int)(2.0 * M_PI * (r + (double)i * 5.0) / 5.0);
		double v = sqrt(a_c * r);
		for (int j = 0; j < j_max; j++) {
			boundary[current + j + NPTS][0] = (r + (double)i * 5.0) * cos(2 * M_PI * j / j_max);
			boundary[current + j + NPTS][1] = (r + (double)i * 5.0) * sin(2 * M_PI * j / j_max);
			boundary[current + j + NPTS][2] = -sin(2 * M_PI * j / j_max)*v;
			boundary[current + j + NPTS][3] = cos(2 * M_PI * j / j_max)*v;
			boundary[current + j + NPTS][4] = 1;
			boundary[current + j + NPTS][5] = 1;
			boundary[current + j + NPTS][6] = 1;
			boundary[current + j + NPTS][7] = 1.0;
		}
		current += j_max;
	}
	r = in_r - 2.5;
	for (int i = 0; i < inside_stages; i++) {
		int j_max = (int)(2.0 * M_PI * (r - (double)i * 5.0) / 5.0);
		double v = -sqrt(1*a_c * r);
		for (int j = 0; j < j_max; j++) {
			boundary[current + j + NPTS][0] = (r - (double)i * 5.0) * cos(2 * M_PI * j / j_max);
			boundary[current + j + NPTS][1] = (r - (double)i * 5.0) * sin(2 * M_PI * j / j_max);
			boundary[current + j + NPTS][2] = -sin(2 * M_PI * j / j_max)*v;
			boundary[current + j + NPTS][3] = cos(2 * M_PI * j / j_max)*v;
			boundary[current + j + NPTS][4] = 1;
			boundary[current + j + NPTS][5] = 1;
			boundary[current + j + NPTS][6] = 1;
			boundary[current + j + NPTS][7] = 1.0;
		}
		current += j_max;
	}
}

void initial_condition(GLfloat(*sup_data)[3], double volume, int total) {
	printf("%f\n", MASS / volume);
	printf("%f\n", MASS / NPTS);
	for (int i = 0; i < NPTS; i++) {
		sup_data[i][0] = MASS/1/volume;
		sup_data[i][1] = MASS/1/NPTS;
		sup_data[i][2] = 293.15;
	}
	for (int i = NPTS; i < (NPTS + total); i++) {
		sup_data[i][0] = MASS/1 / volume;
		sup_data[i][1] = MASS/1 / NPTS;
		sup_data[i][2] = 315;
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

	double radius_min = 20;
	double radius_max = 100.0 - 20;

	double r_max = radius_max + 2.5;
	double a_c = VELOCITY / r_max;
	int total = 0;
	int outside_stages = floor(options->kh / 5);
	for (int i = 0; i < outside_stages; i++)
		total += (int)(2.0 * M_PI * (r_max + (double)i * 5.0) / 5.0);
	double out = total;
	double r_min = radius_min - 2.5;
	int inside_stages = floor(fmin(options->kh, radius_min) / 5);
	for (int i = 0; i < inside_stages; i++)
		total += (int)(2.0 * M_PI * (r_min - (double)i * 5.0) / 5.0);

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
	fillData(data, radius_max-15,radius_min+15);
	boundary_particles_create(outside_stages, inside_stages, data, radius_max, radius_min,a_c);
	initial_condition(sup_data,M_PI*(pow(radius_max,2)- pow(radius_min,2)),total);

	GLfloat(*inner_boundary)[2] = malloc(sizeof(inner_boundary[0]) * 100);
	CHECK_MALLOC(inner_boundary);
	for (int i = 0; i < 100; i++) {
		inner_boundary[i][0] = radius_min * cos(i * 2*M_PI / 100);
		inner_boundary[i][1] = radius_min * sin(i * 2*M_PI / 100);
	}
	GLfloat(*outer_boundary)[2] = malloc(sizeof(outer_boundary[0]) * 300);
	CHECK_MALLOC(outer_boundary);
	for (int i = 0; i < 300; i++) {
		outer_boundary[i][0] = radius_max * cos(i * 2*M_PI / 300);
		outer_boundary[i][1] = radius_max * sin(i * 2*M_PI / 300);
	}

	/* send data to GPU, and receive reference to those data in a points object */
	bov_points_t* particles = bov_particles_new(data, NPTS+total, GL_STATIC_DRAW);
	bov_points_t* inner_particles = bov_points_new(inner_boundary, 100, GL_STATIC_DRAW);
	bov_points_t* outer_particles = bov_points_new(outer_boundary, 300, GL_STATIC_DRAW);

	/* setting particles appearance */
	bov_points_set_width(particles, 0.02);
	bov_points_set_outline_width(particles, 0.0025);
	bov_points_set_width(inner_particles, 0.005);
	bov_points_set_outline_width(inner_particles, 0.0025);
	bov_points_set_width(outer_particles, 0.005);
	bov_points_set_outline_width(outer_particles, 0.0025);

	/* fit particles in a [-0.8 0.8]x[-0.8 0.8] square */
	bov_points_scale(particles, (GLfloat[2]) { 0.008, 0.008 });
	bov_points_set_pos(particles, (GLfloat[2]) { 0.0, -0.1 });
	bov_points_scale(inner_particles, (GLfloat[2]) { 0.008, 0.008 });
	bov_points_set_pos(inner_particles, (GLfloat[2]) { 0.0, -0.1 });
	bov_points_scale(outer_particles, (GLfloat[2]) { 0.008, 0.008 });
	bov_points_set_pos(outer_particles, (GLfloat[2]) { 0.0, -0.1 });
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
			update_predictor(data, sup_data, new_data, nh, timestep, options->kh, grad_w_lucy, 1);
			update_predictor(data, sup_data, new_data, nh, timestep, options->kh, grad_w_lucy, 2);
			timestep /= 2;
			for (int i = 0; i < NPTS; i++) {
				//printf("%f %f %f %f\n", data[i][0], data[i][1], new_data[i][4], new_data[i][5]);
				//sup_data[i][0] = new_data[i][0];
				data[i][2] = new_data[i][2] + new_data[i][4] / sup_data[i][1] * timestep;
				data[i][3] = new_data[i][3] + new_data[i][5] / sup_data[i][1] * timestep;
				if (sqrt(pow(data[i][2], 2) + pow(data[i][3], 2)) > VELOCITY) {
					data[i][2] *= VELOCITY / sqrt(pow(data[i][2], 2) + pow(data[i][3], 2));
					data[i][3] *= VELOCITY / sqrt(pow(data[i][2], 2) + pow(data[i][3], 2));
				}
				data[i][0] += timestep * (new_data[i][2] + new_data[i][4] / sup_data[i][1] * timestep );
				data[i][1] += timestep * (new_data[i][3] + new_data[i][5] / sup_data[i][1] * timestep );
				velocity_to_colormap(sqrt(pow(data[i][2], 2) + pow(data[i][3], 2)), &data[i][4]);
				velocity_to_colormap(sqrt(pow(data[i][2], 2) + pow(data[i][3], 2))* cos(acos(data[i + NPTS][0] / sqrt(pow(data[i + NPTS][0], 2) + pow(data[i + NPTS][1], 2))) - M_PI / 2 + atan(data[i][3] / data[i][2])), &data[i][4]);
				double r = sqrt(pow(data[i][0], 2) + pow(data[i][1], 2));
				double cos_theta = data[i][0] / r;
				double sin_theta = data[i][1] / r;
				double tan_theta_v = data[i][3] / data[i][2];
				double v = sqrt(pow(data[i][2], 2) + pow(data[i][3], 2));
				double cos_theta_v = data[i][2] / v;
				double sin_theta_v = data[i][3] / v;
				if(cos_theta<0 && cos_theta_v > 0)
				velocity_to_colormap(v* sin(atan(sin_theta / cos_theta) - atan(sin_theta_v / cos_theta_v)+M_PI), &data[i][4]);
				else if (cos_theta > 0 && cos_theta_v < 0)
					velocity_to_colormap(v * sin(atan(sin_theta / cos_theta) - atan(sin_theta_v / cos_theta_v) - M_PI), &data[i][4]);
				else
					velocity_to_colormap(v* sin(atan(sin_theta / cos_theta) - atan(sin_theta_v / cos_theta_v)), &data[i][4]);
			}
		}
		neighborhood_update(options, nh, data, 0,total);
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
			double v = -sqrt(1 * a_c * r);
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
		bov_particles_update(particles, data, NPTS+total);
		bov_particles_draw(window, particles, 0, BOV_TILL_END);
		bov_line_loop_draw(window, inner_particles, 0, BOV_TILL_END);
		bov_line_loop_draw(window, outer_particles, 0, BOV_TILL_END);

		bov_text_draw(window, msg);

		// In your actual project, don't wait for events => bov_window_update(window)
		bov_window_update(window);
		prev_t = init_t;
		init_t = clock();
	}
	
	while (!bov_window_should_close(window)) {
		bov_particles_draw(window, particles, 0, BOV_TILL_END);
		bov_line_loop_draw(window, inner_particles, 0, BOV_TILL_END);
		bov_line_loop_draw(window, outer_particles, 0, BOV_TILL_END);
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
