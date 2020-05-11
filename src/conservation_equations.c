#include "conservation_equations.h"

void update_euler(GLfloat(*data)[8], GLfloat(*sup_data)[3], GLfloat(*new_data)[6], neighborhood* nh, double dt, double kh, double (*grad_fun)(double, double, GLfloat**, int, int, int), int type)
{
	if (type == 1) {
		drhodt(data, sup_data, new_data, nh, kh, grad_fun);
	}
	else if (type == 2) {
		dvdt(data, sup_data, new_data, nh, kh, grad_fun);
	}
	else if (type == 3) {
		dTdt(data, sup_data, new_data, nh, kh, grad_fun);
	}
	for (int i = 0; i < NPTS; i++) {
		if (type == 1) {
			//new_data[nh[i].index][0] = sup_data[nh[i].index][0] + new_data[nh[i].index][0] * dt;
		}
		else if (type == 2) {
			//printf("%f %f\n", new_data[nh[i].index][2], new_data[nh[i].index][3]);
			new_data[nh[i].index][2] = data[nh[i].index][2] + new_data[nh[i].index][2] * dt;
			new_data[nh[i].index][3] = data[nh[i].index][3] + new_data[nh[i].index][3] * dt;
		}
		else if (type==3)
			new_data[nh[i].index][1] = sup_data[nh[i].index][2] + new_data[nh[i].index][1] * dt;
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
	else if (type == 3) {
		dTdt(data, sup_data, new_data, nh, kh, grad_fun);
	}
	for (int i = 0; i < NPTS; i++) {
		if (type == 1) {
			new_data[nh[i].index][0] = sup_data[nh[i].index][0] + new_data[nh[i].index][0] * dt;
		}
		else if (type == 2) {
			//printf("%f %f\n", new_data[nh[i].index][2], new_data[nh[i].index][3]);
			new_data[nh[i].index][2] = data[nh[i].index][2] + new_data[nh[i].index][2] * dt/2;
			new_data[nh[i].index][3] = data[nh[i].index][3] + new_data[nh[i].index][3] * dt/2;
			data[i][2] = new_data[i][2] + new_data[i][4] / sup_data[i][1] * dt/2;
			data[i][3] = new_data[i][3] + new_data[i][5] / sup_data[i][1] * dt / 2;
			data[i][0] += dt/2 * data[i][2];
			data[i][1] += dt/2 * data[i][3];
		}
		else if (type == 3)
			sup_data[nh[i].index][2] = sup_data[nh[i].index][2] + new_data[nh[i].index][1] * dt/2;
	}
	if (type == 1) {
		drhodt(data, sup_data, new_data, nh, kh, grad_fun);
	}
	else if (type == 2) {
		dvdt(data, sup_data, new_data, nh, kh, grad_fun);
	}
	else if (type == 3) {
		dTdt(data, sup_data, new_data, nh, kh, grad_fun);
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
		if (type == 3) {
			new_data[nh[i].index][1] = sup_data[nh[i].index][2] + new_data[nh[i].index][1] * dt/2;
		}
	}
}

double dyn_pressure(double density) {
	return 1100 * INITDENSITY / 7 * (pow(density / INITDENSITY, 7) - 1);
}

void drhodt(GLfloat(*data)[8], GLfloat(*sup_data)[3], GLfloat(*new_data)[6], neighborhood* nh, double kh, double (*grad_fun)(double, double, GLfloat**, int, int, int)) {
	for (int i = 0; i < NPTS; i++) {
		new_data[nh[i].index][0] = 0;
		double den = 0;
		for (neighbours* j = nh[nh[i].index].list; j; j = j->next) {
			//printf("i %i j %i v %f v %f g %f g %f\n",i, j->index, data[nh[i].index][2], data[j->index][2], grad_fun(j->distance, kh, data, nh[i].index, j->index, 1), data[nh[i].index][3], data[j->index][3], grad_fun(j->distance, kh, data, nh[i].index, j->index, 0));
			//new_data[nh[i].index][0] += sup_data[j->index][1] * ((data[nh[i].index][2] - data[j->index][2]) * grad_fun(j->distance, kh, data, nh[i].index, j->index, 1) + (data[nh[i].index][3] - data[j->index][3]) * grad_fun(j->distance, kh, data, nh[i].index, j->index, 0));
			new_data[nh[i].index][0] += sup_data[j->index][1] * grad_fun(j->distance, kh, data, nh[i].index, j->index, 1);
			den += grad_fun(j->distance, kh, data, nh[i].index, j->index, 1) * sup_data[j->index][1] / sup_data[j->index][0];
		}
		new_data[nh[i].index][0] /= den;
		//printf("%f\n", new_data[nh[i].index][0]);
	}
}

void dvdt(GLfloat(*data)[8], GLfloat(*sup_data)[3], GLfloat(*new_data)[6], neighborhood* nh, double kh, double (*grad_fun)(double, double, GLfloat**, int, int, int)) {
	double* dPdx = calloc(NPTS, sizeof(double));
	double* dPdy = calloc(NPTS,sizeof(double));
	for (int i = 0; i < NPTS; i++) {
		double denx = 0;
		double deny = 0;
		for (neighbours* j = nh[nh[i].index].list; j; j = j->next) {
			dPdx[nh[i].index] += (dyn_pressure(sup_data[nh[i].index][0]) / pow(sup_data[nh[i].index][0], 2) + dyn_pressure(sup_data[j->index][0]) / pow(sup_data[j->index][0], 2)) * grad_fun(j->distance, kh, data, nh[i].index, j->index, 1);
			dPdy[nh[i].index] += (dyn_pressure(sup_data[nh[i].index][0]) / pow(sup_data[nh[i].index][0], 2) + dyn_pressure(sup_data[j->index][0]) / pow(sup_data[j->index][0], 2)) * grad_fun(j->distance, kh, data, nh[i].index, j->index, 0);
			denx += grad_fun(j->distance, kh, data, nh[i].index, j->index, 1)* (data[j->index][0]-data[nh[i].index][0])* sup_data[j->index][1]/ sup_data[j->index][0];
			deny += grad_fun(j->distance, kh, data, nh[i].index, j->index, 0) * (data[j->index][1] - data[nh[i].index][1]) * sup_data[j->index][1] / sup_data[j->index][0];
		}
		dPdx[nh[i].index] /= denx;
		dPdy[nh[i].index] /= deny;
	}
	for (int i = 0; i < NPTS; i++) {
		int gx = 0;
		int gy = -0;
		new_data[nh[i].index][2] = gx;
		new_data[nh[i].index][3] = gy;
		new_data[nh[i].index][4] = 0;
		new_data[nh[i].index][5] = 0;
		for (neighbours* j = nh[nh[i].index].list; j; j = j->next) {
			if (j->index >= NPTS && j->distance<=R_0) {
				new_data[nh[i].index][4] += PSY / pow(j->distance, 2) * (pow(R_0 / j->distance, 12) - pow(R_0 / j->distance, 4)) * (data[nh[i].index][0] - data[j->index][0]);
				new_data[nh[i].index][5] += PSY / pow(j->distance, 2) * (pow(R_0 / j->distance, 12) - pow(R_0 / j->distance, 4)) * (data[nh[i].index][1] - data[j->index][1]);
			}
			new_data[nh[i].index][2] += sup_data[j->index][1] * dPdx[nh[i].index] + 2 * MU / sup_data[nh[i].index][0] * sup_data[j->index][1] / sup_data[j->index][0] * (data[nh[i].index][2] - data[j->index][2]) * ((data[nh[i].index][0] - data[j->index][0]) * grad_fun(j->distance, kh, data, nh[i].index, j->index, 1) + (data[nh[i].index][1] - data[j->index][1]) * grad_fun(j->distance, kh, data, nh[i].index, j->index, 0)) / pow(j->distance, 2);
			new_data[nh[i].index][3] += sup_data[j->index][1] * dPdy[nh[i].index] + 2 * MU / sup_data[nh[i].index][0] * sup_data[j->index][1] / sup_data[j->index][0] * (data[nh[i].index][3] - data[j->index][3]) * ((data[nh[i].index][0] - data[j->index][0]) * grad_fun(j->distance, kh, data, nh[i].index, j->index, 1) + (data[nh[i].index][1] - data[j->index][1]) * grad_fun(j->distance, kh, data, nh[i].index, j->index, 0)) / pow(j->distance, 2);
		}
	}
	free(dPdx);
	free(dPdy);
}

void dTdt(GLfloat(*data)[8], GLfloat(*sup_data)[3], GLfloat(*new_data)[6], neighborhood* nh, double kh, double (*grad_fun)(double, double, GLfloat**, int, int, int)) {
	for (int i = 0; i < NPTS; i++) {
		new_data[nh[i].index][1] = 0;
		for (neighbours* j = nh[nh[i].index].list; j; j = j->next) {
			//printf("i %i j %i v %f v %f g %f g %f\n",i, j->index, data[nh[i].index][2], data[j->index][2], grad_fun(j->distance, kh, data, nh[i].index, j->index, 1), data[nh[i].index][3], data[j->index][3], grad_fun(j->distance, kh, data, nh[i].index, j->index, 0));
			new_data[nh[i].index][1] += -dyn_pressure(sup_data[nh[i].index][0])*sup_data[j->index][1] * ((data[nh[i].index][2] - data[j->index][2]) * grad_fun(j->distance, kh, data, nh[i].index, j->index, 1) + (data[nh[i].index][3] - data[j->index][3]) * grad_fun(j->distance, kh, data, nh[i].index, j->index, 0))+alph* 2 * sup_data[j->index][1] / sup_data[j->index][0] * (sup_data[nh[i].index][2] - sup_data[j->index][2]) * ((data[nh[i].index][0] - data[j->index][0]) * grad_fun(j->distance, kh, data, nh[i].index, j->index, 1) + (data[nh[i].index][1] - data[j->index][1]) * grad_fun(j->distance, kh, data, nh[i].index, j->index, 0)) / pow(j->distance, 2);
		}
	}
}