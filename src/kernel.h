#ifndef KERNEL_H
#define KERNEL_H


#include "BOV.h"

#include "neighborhood_search.h"
#include "shared_variables.h"
#include <time.h>
#include <math.h>

typedef struct neighborhood neighborhood;

/*
 Implementation of the gradient of the kernel cubic spline function
 Input : the distance between the particles, the radius of the neighborhood and the distance between particle in x or y direction regarding the desired weight
 Output : the coefficient that represents the importance of each particles on the other particles.
 */
double grad_w_cubic(double distance, double kh, GLfloat(*data)[8], int i, int j, int is_x);

/*
 Implementation of the gradient of the kernel quintic spline function
 Input : the distance between the particles, the radius of the neighborhood and the distance between particle in x or y direction regarding the desired weight
 Output : the coefficient that represents the importance of each particles on the other particles.
 */
double grad_w_quinticspline(double distance, double kh, GLfloat(*data)[8], int i, int j, int is_x);

/*
 Implementation of the gradient of the kernel new quartic spline function
 Input : the distance between the particles, the radius of the neighborhood and the distance between particle in x or y direction regarding the desired weight
 Output : the coefficient that represents the importance of each particles on the other particles.
 */
double grad_w_newquartic(double distance, double kh, GLfloat(*data)[8], int i, int j, int is_x);

/*
 Implementation of the gradient of the Lucy quartic kernel function
 Input : the distance between the particles, the radius of the neighborhood and the distance between particle in x or y direction regarding the desired weight
 Output : the coefficient that represents the importance of each particles on the other particles.
 */
double grad_w_lucy(double distance, double kh, GLfloat(*data)[8], int i, int j, int is_x);

double w_cubic(double distance, double kh);

double w_quinticspline(double distance, double kh);

double w_newquartic(double distance, double kh);

double w_lucy(double distance, double kh, GLfloat(*data)[8], int i, int j, int is_x);

#endif
