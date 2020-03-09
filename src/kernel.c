#include "kernel.h"
#include <math.h>


// Implementation of the kernel cubic function and return the weight regarding the distance and the radius of the circle


void kernel(GLfloat(*data)[14], neighborhood* nh, double kh) {
    for (int i = 0; i < NPTS; i++) {
        double val_node_x = data[i][8];
        double val_node_y = data[i][9];
        int nNeigh = nh[i].nNeighbours;
        double val_div = 0;
        double val_grad_x = 0;
        double val_grad_y = 0;
        double val_lapl = 0;
        double dens2 = pow(DENSITY, 2);
        neighbours* List = nh[i].list;
        if (nNeigh > 0) {
            for (int j = 0; j < nNeigh; j++) {
                int index_node2 = List->index;
                double distance = List->distance;
                double d_x = data[index_node2][0] - data[i][0];
                double d_y = data[index_node2][1] - data[i][1];
                
                /*
                 You can choose here the desired kernel function for your code.
                 */
                
                //double weight_x = grad_w_cubic(distance, kh, d_x);
                //double weight_y = grad_w_cubic(distance, kh, d_y);
                
                double weight_x = grad_w_lucy(distance, kh, d_x);
                double weight_y = grad_w_lucy(distance, kh, d_y);
                
                //double weight_x = grad_w_newquartic(distance, kh, d_x);
                //double weight_y = grad_w_newquartic(distance, kh, d_y);
                
                //double weight_x = grad_w_quinticspline(distance, kh, d_x);
                //double weight_y = grad_w_quinticspline(distance, kh, d_y);
                
                val_div += -MASS / DENSITY * ((data[index_node2][8] - val_node_x) * weight_x + (data[index_node2][9] - val_node_y) * weight_y);
                val_grad_x += -DENSITY * MASS * ((val_node_x / dens2) + (data[index_node2][8] / dens2)) * weight_x;
                val_grad_y += -DENSITY * MASS * ((val_node_x / dens2) + (data[index_node2][8] / dens2)) * weight_y;
                val_lapl += 2.0 * MASS / DENSITY * (val_node_x - data[index_node2][8]) * (d_x * weight_x + d_y * weight_y) / pow(distance,2);
                
                List = List->next;
            }
        }
        // All the values of the divergent gradient and laplacien are stored in the data table
        data[i][10] = val_div;
        data[i][11] = val_grad_x;
        data[i][12] = val_grad_y;
        data[i][13] = val_lapl;
    }
    
    //Computation of the error based on the already know function.
    for (int j = 0; j < NPTS; j++) {
        double exact = 3 * pow(data[j][0], 2);
        double error = exact - data[j][10];
    }
}

/*
 Implementation of the kernel cubic function and return the weight regarding the distance and the radius of the circle
 */
double grad_w_cubic(double distance, double kh, double d)
{
    double h = kh / 2;
    double q = distance / h;
    double weight = 0;
    double alpha_d = 15 / (7 * M_PI * pow(h, 2));
    
    if (q > 0) {
        if (q <= 1) {
            weight = ((-2.0 * q + 1.5 * pow(q, 2)) * d) / (pow(h, 2) * q );
        }
        else {
            if (q > 2) {
                weight = 0;
                printf("too big distance ");
            }
            else {
                weight = (-0.5 * pow((2 - q), 2) * d) / (pow(h, 2) * q);
            }
        }
    }
    return weight*alpha_d;
}

double grad_w_lucy(double distance, double kh, double d)
{
    double h = kh / 1;
    double q = distance / h;
    double grad_w = 0;
    double alpha_d = (5 / (M_PI * pow(h, 2)));
    if (q >= 0 && q <= 1)
    {
        grad_w = (-12.0 * q * alpha_d * pow((1 - q), 2)) * d /(pow(h,2) * q);
    }
    else
    {
        grad_w = 0.0;
    }
    return  grad_w;
}


double grad_w_newquartic(double distance, double kh, double d)
{
    double h = kh / 2;
    double q = distance / h;
    double grad_w = 0;
    double alpha_d = (15 / (7 * M_PI * pow(h, 2)));
    if (q >= 0 && q <= 2)
    {
        grad_w = (alpha_d * (-(9.0 / 4.0) * q + (19.0 / 8.0) * pow(q, 2) - (5.0 / 8.0) * pow(q, 3)) * d) / (pow(h, 2) * q);
    }
    else
    {
        grad_w = 0;
    }
    return  grad_w;
}

double grad_w_quinticspline(double distance, double kh, double d)
{
    double h = kh / 3;
    double q = distance / h;
    double x_x, y_y;
    double grad_w = 0;
    double alpha_d = (7 / (478 * M_PI * pow(h, 2)));
    double dq = d / (h * distance);
    if (q >= 0 && q <= 1)
    {
        grad_w = (alpha_d * (-5 * pow((3 - q), 4) + 30 * pow((2 - q), 4) - 75 * pow((1 - q), 4)) * d) / (pow(h, 2) * q);
    }
    else  if (q > 1 && q <= 2)
    {
        grad_w = (alpha_d * (-5 * pow((3 - q), 4) + 30 * pow((2 - q), 4)) * d) / (pow(h, 2) * q);
    }
    else  if (q > 2 && q <= 3)
    {
        grad_w = (alpha_d * (-5 * pow((3 - q), 4)) * d) / (pow(h, 2) * q);
    }
    else
    {
        grad_w = 0;
    }
    return  grad_w;
}




