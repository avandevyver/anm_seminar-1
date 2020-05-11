#ifndef CONSERVATION_EQUATIONS_H
#define CONSERVATION_EQUATIONS_H

#include "BOV.h"
#include <time.h>
#include <math.h>
#include "shared_variables.h"
#include "neighborhood_search.h"

double dyn_pressure(double density);

void update_euler(GLfloat(*data)[8], GLfloat(*sup_data)[3], GLfloat(*new_data)[6], neighborhood* nh, double dt, double kh, double (*grad_fun)(double, double, GLfloat**, int, int, int), int type);

void update_predictor(GLfloat(*data)[8], GLfloat(*sup_data)[3], GLfloat(*new_data)[6], neighborhood* nh, double dt, double kh, double (*grad_fun)(double, double, GLfloat**, int, int, int), int type);

void drhodt(GLfloat(*data)[8], GLfloat(*sup_data)[3], GLfloat(*new_data)[6], neighborhood* nh, double kh, double (*grad_fun)(double, double, GLfloat**, int, int, int));

void dvdt(GLfloat(*data)[8], GLfloat(*sup_data)[3], GLfloat(*new_data)[6], neighborhood* nh, double kh, double (*grad_fun)(double, double, GLfloat**, int, int, int));

void dTdt(GLfloat(*data)[8], GLfloat(*sup_data)[3], GLfloat(*new_data)[6], neighborhood* nh, double kh, double (*grad_fun)(double, double, GLfloat**, int, int, int));
#endif