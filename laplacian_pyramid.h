#ifndef LAPLACIAN_PYRAMID_H
#define LAPLACIAN_PYRAMID_H
using namespace std;
#include <iostream>
#include <stdio.h>
#include <fstream>
#include <math.h>
#include <string.h>
#include <stdlib.h>

double **laplacian_pyramid(double *I, double *filter, int *dims, int nlev, int *subwindow, int **dimsout);//return subwindow_child too
#endif 