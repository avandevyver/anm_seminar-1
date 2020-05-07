#ifndef CONSERVATION_EQUATIONS_H
#define CONSERVATION_EQUATIONS_H

#include "BOV.h"
#include <time.h>
#include <math.h>
#include "neighborhood_search.h"

#define INIT_DENSITY 1000
#define INIT_TEMPERATURE 293.15
#define BOUNDARY_VELOCITY 1

void update_euler(GLfloat(*data)[8], GLfloat(*sup_data)[3], GLfloat(*new_data)[6], neighborhood* nh, double dt, double kh, double (*grad_fun)(double, double, GLfloat**, int, int, int), int type);

void initialise_leapfrog(GLfloat(*data)[8], GLfloat(*coord)[2], float(*sol)[5], double dt, double kh, neighborhood* nhList, double sourceTemp, float* newSolm);

void update_leapfrog(GLfloat(*data)[8], GLfloat(*coord)[2], float(*sol)[5], double dt, double kh, neighborhood* nhList, double sourceTemp, float* newSolm);

void update_predictor(GLfloat(*data)[8], GLfloat(*sup_data)[3], GLfloat(*new_data)[6], neighborhood* nh, double dt, double kh, double (*grad_fun)(double, double, GLfloat**, int, int, int), int type);

void update_simple(GLfloat(*data)[8], GLfloat(*coord)[2], float(*sol)[5], float dt, double kh, neighborhood* nhList, double sourceTemp);

void fillDataSem3(GLfloat(*data)[8], GLfloat(*coord)[2], int nPoints, double minTemp, double maxTemp, float(*sol)[5]);

void fillDataSem3_zero(GLfloat(*data)[8], GLfloat(*coord)[2], int nPoints, double minTemp, double maxTemp, float(*sol)[5]);

void drhodt(GLfloat(*data)[8], GLfloat(*sup_data)[3], GLfloat(*new_data)[6], neighborhood* nh, double kh, double (*grad_fun)(double, double, GLfloat**, int, int, int));

void dvdt(GLfloat(*data)[8], GLfloat(*sup_data)[3], GLfloat(*new_data)[6], neighborhood* nh, double kh, double (*grad_fun)(double, double, GLfloat**, int, int, int));

void tempToColor(float temp, float color[3], double initTemp, double sourceTemp);
#endif