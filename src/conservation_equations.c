#include "conservation_equations.h"

#define INIT_DENSITY 0.000027
#define MU 0.015
#define PSY 3
#define R_0 10
/*
Update the temperature of all points (value of temperature of point i is in sol[i][0])
*/
void update_euler(GLfloat(*data)[8], GLfloat(*sup_data)[3], GLfloat(*new_data)[6], neighborhood* nh, double dt, double kh, double (*grad_fun)(double, double, GLfloat**, int, int, int), int type)
{
	if (type == 1) {
		drhodt(data, sup_data, new_data, nh, kh, grad_fun);
	}
	else if (type == 2) {
		dvdt(data, sup_data, new_data, nh, kh, grad_fun);
	}
	for (int i = 0; i < NPTS; i++) {
		if (type == 1) {
			new_data[nh[i].index][0] = sup_data[nh[i].index][0] + new_data[nh[i].index][0] * dt;
		}
		else if (type == 2) {
			//printf("%f %f\n", new_data[nh[i].index][2], new_data[nh[i].index][3]);
			new_data[nh[i].index][2] = data[nh[i].index][2] + new_data[nh[i].index][2] * dt;
			new_data[nh[i].index][3] = data[nh[i].index][3] + new_data[nh[i].index][3] * dt;
		}
	}
	/*for (int i = 0; i < NPTS; i++) {
		if (newSol[i] < 374 && newSol[i]>272) {
			sol[i][0] = newSol[i]; // condition permettant d'éviter certaines instabilités
		}
	}*/
}

void update_predictor(GLfloat(*data)[8], GLfloat(*sup_data)[3], GLfloat(*new_data)[6], neighborhood* nh, double dt, double kh, double (*grad_fun)(double, double, GLfloat**, int, int, int), int type)
{
	if (type == 1) {
		drhodt(data, sup_data, new_data, nh, kh, grad_fun);
	}
	else if (type == 2) {
		dvdt(data, sup_data, new_data, nh, kh, grad_fun);
	}
	for (int i = 0; i < NPTS; i++) {
		if (type == 1) {
			new_data[nh[i].index][0] = sup_data[nh[i].index][0] + new_data[nh[i].index][0] * dt;
		}
		else if (type == 2) {
			//printf("%f %f\n", new_data[nh[i].index][2], new_data[nh[i].index][3]);
			new_data[nh[i].index][2] = data[nh[i].index][2] + new_data[nh[i].index][2] * dt/2;
			new_data[nh[i].index][3] = data[nh[i].index][3] + new_data[nh[i].index][3] * dt/2;
			data[i][2] = new_data[i][2] + new_data[i][4] / sup_data[i][1] * dt / 2;
			data[i][3] = new_data[i][3] + new_data[i][5] / sup_data[i][1] * dt / 2;
			data[i][0] += dt/2 * data[i][2];
			data[i][1] += dt/2 * data[i][3];
		}
	}
	if (type == 1) {
		drhodt(data, sup_data, new_data, nh, kh, grad_fun);
	}
	else if (type == 2) {
		dvdt(data, sup_data, new_data, nh, kh, grad_fun);
	}
	for (int i = 0; i < NPTS; i++) {
		if (type == 1) {
			new_data[nh[i].index][0] = sup_data[nh[i].index][0] + new_data[nh[i].index][0] * dt;
		}
		else if (type == 2) {
			//printf("%f %f\n", new_data[nh[i].index][2], new_data[nh[i].index][3]);
			new_data[nh[i].index][2] = data[nh[i].index][2] + new_data[nh[i].index][2] * dt / 2;
			new_data[nh[i].index][3] = data[nh[i].index][3] + new_data[nh[i].index][3] * dt / 2;
		}
	}
}

// Leapfrog
void initialise_leapfrog(GLfloat(*data)[8], GLfloat(*coord)[2], float(*sol)[5], double dt, double kh, neighborhood* nhList, double sourceTemp, float newSolm[1000])
{
	float newSol[1000];
	for (int i = 0; i < NPTS; i++) {
		newSol[i] = sol[i][0];
		if (sol[i][0] != sourceTemp) {
			neighbours* point = nhList[i].list;
			for (int j = 0; j < (int)nhList[i].nNeighbours; j++) {
				float q = point->distance / kh;
				float gradW = 5 * (-12 * q + 24 * q * q - 12 * q * q * q) / point->distance;
				//float gradW = 1/(point->distance
				int numbers = (int)nhList[point->index].nNeighbours;
				float to_add = 2 * (sol[i][0] - sol[point->index][0]) * gradW / (point->distance * numbers);
				newSol[i] += -0.5 * dt * to_add;
				point = point->next;
			}
		}
	}
	for (int i = 0; i < NPTS; i++) {
		newSolm[i] = newSol[i];
	}
}

void update_leapfrog(GLfloat(*data)[8], GLfloat(*coord)[2], float(*sol)[5], double dt, double kh, neighborhood* nhList, double sourceTemp, float newSolm[1000])
{
	float newSolp[1000];
	for (int i = 0; i < NPTS; i++) {
		newSolp[i] = newSolm[i];
		if (sol[i][0] != sourceTemp) {
			neighbours* point = nhList[i].list;
			for (int j = 0; j < (int)nhList[i].nNeighbours; j++) {
				float q = point->distance / kh;
				float gradW = 5 * (-12 * q + 24 * q * q - 12 * q * q * q) / point->distance;
				//float gradW = 1/(point->distance
				int numbers = (int)nhList[point->index].nNeighbours;
				float to_add = 2 * (newSolm[i] - newSolm[point->index]) * gradW / (point->distance * numbers);
				newSolp[i] += dt * to_add;
				point = point->next;
			}
		}
	}
	for (int i = 0; i < NPTS; i++) {
		if ((0.5 * (newSolp[i] + newSolm[i])) < 374 && 0.5 * (newSolp[i] + newSolm[i]) > 272) {
			sol[i][0] = 0.5 * (newSolp[i] + newSolm[i]);
			newSolm[i] = newSolp[i];
		}

	}
}

// PREDICTOR CORRECTOR

/*void update_predictor(GLfloat(*data)[8], GLfloat(*coord)[2], float(*sol)[5], double dt, double kh, neighborhood* nhList, double sourceTemp)
{
	float newSol[1000];
	float laplacian1[1000];
	float laplacian2[1000];
	for (int i = 0; i < NPTS; i++) {
		newSol[i] = sol[i][0];
		laplacian1[i] = 0;
		if (sol[i][0] != sourceTemp) {
			neighbours* point = nhList[i].list;
			for (int j = 0; j < (int)nhList[i].nNeighbours; j++) {
				float q = point->distance / kh;
				float gradW = 5 * (-12 * q + 24 * q * q - 12 * q * q * q) / point->distance;
				//float gradW = 1/(point->distance
				int numbers = (int)nhList[point->index].nNeighbours;
				float to_add = 2 * (sol[i][0] - sol[point->index][0]) * gradW / (point->distance * numbers);
				newSol[i] += 0.5 * dt * to_add;
				laplacian1[i] += dt * to_add;
				point = point->next;
			}
		}
	}
	for (int i = 0; i < NPTS; i++) {
		laplacian2[i] = 0;
		if (sol[i][0] != sourceTemp) {
			neighbours* point = nhList[i].list;
			for (int j = 0; j < (int)nhList[i].nNeighbours; j++) {
				float q = point->distance / kh;
				float gradW = 5 * (-12 * q + 24 * q * q - 12 * q * q * q) / point->distance;
				//float gradW = 1/(point->distance
				int numbers = (int)nhList[point->index].nNeighbours;
				float to_add = 2 * (newSol[i] - newSol[point->index]) * gradW / (point->distance * numbers);
				laplacian2[i] += dt * to_add;
				point = point->next;
			}
		}
	}
	for (int i = 0; i < NPTS; i++) {
		if ((sol[i][0] + 0.5 * (laplacian1[i] + laplacian2[i])) < 374 && (sol[i][0] + 0.5 * (laplacian1[i] + laplacian2[i])) > 272) {
			sol[i][0] += 0.5 * (laplacian1[i] + laplacian2[i]); // renforce la stabilité
		}
	}
}*/



/*
Update the temperature of all points (value of temperature of point i is in sol[i][0])
*/
void update_simple(GLfloat(*data)[8], GLfloat(*coord)[2], float(*sol)[5], float dt, double kh, neighborhood* nhList, double sourceTemp)
{
	float newSol[1000];
	for (int i = 0; i < NPTS; i++) {
		newSol[i] = sol[i][0];
		if (sol[i][0] != sourceTemp) {
			neighbours* point = nhList[i].list;
			float sumDist = 1;
			for (int j = 0; j < (int)nhList[i].nNeighbours; j++) {
				sumDist += (float)(((kh - point->distance)) * ((kh - point->distance)));
				newSol[i] += (float)((kh - point->distance) * (kh - point->distance) * (double)sol[point->index][0]);
				point = point->next;
			}
			newSol[i] = newSol[i] / sumDist;
		}

	}
	for (int i = 0; i < NPTS; i++) {
		sol[i][0] = newSol[i];
	}
}

void fillDataSem3(GLfloat(*data)[8], GLfloat(*coord)[2], int nPoints, double minTemp, double maxTemp, float(*sol)[5])
{
	for (int i = 0; i < NPTS; i++) {
		//data[i][0] = (i / (int)sqrt(NPTS)) * 200 / sqrt(NPTS) - 100;
		//data[i][1] = (i % (int)sqrt(NPTS)) * 200 / sqrt(NPTS) - 100;
		data[i][0] = rand() * 200.0 / RAND_MAX - 100.0; // x (rand between -100 and 100)
		data[i][1] = rand() * 200.0 / RAND_MAX - 100.0; // y (rand between -100 and 100)
		coord[i][0] = data[i][0];
		coord[i][1] = data[i][1];
		data[i][2] = 0; //fixed particles for Seminar 3
		data[i][3] = 0; //fixed particles for Seminar 3
		//data[i][2] = rand() * 2.0 / RAND_MAX - 1.0; //Random starting speed
		//data[i][3] = rand() * 2.0 / RAND_MAX - 1.0; //Random starting speed
		if (data[i][0] < -90) {
			//tempToColor(maxTemp, &data[i][4], minTemp, maxTemp); // fill color
			sol[i][0] = maxTemp;
		}
		else {
			//tempToColor(minTemp, &data[i][4], minTemp, maxTemp); // fill color
			sol[i][0] = minTemp;
		}
		data[i][7] = 1; // transparency

	}
}

void fillDataSem3_zero(GLfloat(*data)[8], GLfloat(*coord)[2], int nPoints, double minTemp, double maxTemp, float(*sol)[5])
{
	for (int i = 0; i < NPTS; i++) {
		if (data[i][0] < -90) {
			//tempToColor(maxTemp, &data[i][4], minTemp, maxTemp); // fill color
			sol[i][0] = maxTemp;
		}
		else {
			//tempToColor(minTemp, &data[i][4], minTemp, maxTemp); // fill color
			sol[i][0] = minTemp;
		}
		data[i][7] = 1; // transparency
	}
}

double dyn_pressure(double density, double temperature) {
	return 1.4 * 8.314 * temperature * INIT_DENSITY / 7 * (pow(density / INIT_DENSITY, 7) - 1);
	//return 1.4 * 8.314 * temperature * density;
}

void drhodt(GLfloat(*data)[8], GLfloat(*sup_data)[3], GLfloat(*new_data)[6], neighborhood* nh, double kh, double (*grad_fun)(double, double, GLfloat**, int, int, int)) {
	for (int i = 0; i < NPTS; i++) {
		new_data[nh[i].index][0] = 0;
		for (neighbours* j = nh[nh[i].index].list; j; j = j->next) {
			//printf("i %i j %i v %f v %f g %f g %f\n",i, j->index, data[nh[i].index][2], data[j->index][2], grad_fun(j->distance, kh, data, nh[i].index, j->index, 1), data[nh[i].index][3], data[j->index][3], grad_fun(j->distance, kh, data, nh[i].index, j->index, 0));
			new_data[nh[i].index][0] += sup_data[j->index][1] * ((data[nh[i].index][2] - data[j->index][2]) * grad_fun(j->distance, kh, data, nh[i].index, j->index, 1) + (data[nh[i].index][3] - data[j->index][3]) * grad_fun(j->distance, kh, data, nh[i].index, j->index, 0));
		}
	}
}

void dvdt(GLfloat(*data)[8], GLfloat(*sup_data)[3], GLfloat(*new_data)[6], neighborhood* nh, double kh, double (*grad_fun)(double, double, GLfloat**, int, int, int)) {
	for (int i = 0; i < NPTS; i++) {
		int gx = 0;
		int gy = -0;
		new_data[nh[i].index][2] = gx;
		new_data[nh[i].index][3] = gy;
		new_data[nh[i].index][4] = 0;
		new_data[nh[i].index][5] = 0;
		for (neighbours* j = nh[nh[i].index].list; j; j = j->next) {
			//printf("i %i j %i nd %f nd %f m %f rho %f p %f p %f g %f g %f\n",i, j->index, new_data[nh[i].index][2], new_data[nh[i].index][3], sup_data[j->index][1], sup_data[j->index][0], dyn_pressure(sup_data[nh[i].index][0], sup_data[nh[i].index][2]), dyn_pressure(sup_data[j->index][0], sup_data[j->index][2]), grad_fun(j->distance, kh, data, nh[i].index, j->index, 1), grad_fun(j->distance, kh, data, nh[i].index, j->index, 0));
			if (j->index >= NPTS && j->distance<=R_0) {
				new_data[nh[i].index][4] += PSY / pow(j->distance, 2) * (pow(R_0 / j->distance, 12) - pow(R_0 / j->distance, 4)) * (data[nh[i].index][0] - data[j->index][0]);
				new_data[nh[i].index][5] += PSY / pow(j->distance, 2) * (pow(R_0 / j->distance, 12) - pow(R_0 / j->distance, 4)) * (data[nh[i].index][1] - data[j->index][1]);
			}
			new_data[nh[i].index][2] += sup_data[j->index][1] * ((dyn_pressure(sup_data[nh[i].index][0], sup_data[nh[i].index][2]) / pow(sup_data[nh[i].index][0], 2) + dyn_pressure(sup_data[j->index][0], sup_data[j->index][2]) / pow(sup_data[j->index][0], 2)) * grad_fun(j->distance, kh, data, nh[i].index, j->index, 1) + 2 * MU / sup_data[nh[i].index][0] * sup_data[j->index][1] / sup_data[j->index][0] * (data[nh[i].index][2] - data[j->index][2]) * ((data[nh[i].index][0] - data[j->index][0]) * grad_fun(j->distance, kh, data, nh[i].index, j->index, 1) + (data[nh[i].index][1] - data[j->index][1]) * grad_fun(j->distance, kh, data, nh[i].index, j->index, 0)) / pow(j->distance, 2));
			new_data[nh[i].index][3] += sup_data[j->index][1] * ((dyn_pressure(sup_data[nh[i].index][0], sup_data[nh[i].index][2]) / pow(sup_data[nh[i].index][0], 2) + dyn_pressure(sup_data[j->index][0], sup_data[j->index][2]) / pow(sup_data[j->index][0], 2)) * grad_fun(j->distance, kh, data, nh[i].index, j->index, 0) + 2 * MU / sup_data[nh[i].index][0] * sup_data[j->index][1] / sup_data[j->index][0] * (data[nh[i].index][3] - data[j->index][3]) * ((data[nh[i].index][0] - data[j->index][0]) * grad_fun(j->distance, kh, data, nh[i].index, j->index, 1) + (data[nh[i].index][1] - data[j->index][1]) * grad_fun(j->distance, kh, data, nh[i].index, j->index, 0)) / pow(j->distance, 2));
			//new_data[nh[i].index][2] += sup_data[j->index][1] * ( 2 * MU / sup_data[nh[i].index][0] * sup_data[j->index][1] / sup_data[j->index][0] * (data[nh[i].index][2] - data[j->index][2]) * ((data[nh[i].index][0] - data[j->index][0]) * grad_fun(j->distance, kh, data, nh[i].index, j->index, 1) + (data[nh[i].index][1] - data[j->index][1]) * grad_fun(j->distance, kh, data, nh[i].index, j->index, 0)) / pow(j->distance, 2));
			//new_data[nh[i].index][3] += sup_data[j->index][1] * (2 * MU / sup_data[nh[i].index][0] * sup_data[j->index][1] / sup_data[j->index][0] * (data[nh[i].index][3] - data[j->index][3]) * ((data[nh[i].index][0] - data[j->index][0]) * grad_fun(j->distance, kh, data, nh[i].index, j->index, 1) + (data[nh[i].index][1] - data[j->index][1]) * grad_fun(j->distance, kh, data, nh[i].index, j->index, 0)) / pow(j->distance, 2));
		}
	}
}

/*void tempToColor(float temp, float color[3], double initTemp, double sourceTemp)
{
	float v = (temp - initTemp) / (sourceTemp - initTemp);
	/*
	float v1 = 3.5 * (v - 0.7);
	float v2 = 1.25 * v;
	float v3 = fminf(0.5, v) * 2.0;

	color[0] = -v1 * v1 + 1.0f;
	color[1] = 6.0f * v2 * v2 * (1.0f - v2);
	color[2] = 5.5f * v3 * (1.0f - v3) * (1.0f - v3);
	
	// alternative: classical jet colormap
	color[0] = 1.5 - 4.0 * fabs(v - 0.75);
	color[1] = 1.5 - 4.0 * fabs(v - 0.5);
	color[2] = 1.5 - 4.0 * fabs(v - 0.25);
}*/