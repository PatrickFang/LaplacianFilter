#ifndef RECONSTRUCT_LAPLACIAN_PYRAMID_H
#define RECONSTRULAPLACIAN_PYRAMID_H
using namespace std;
#include <iostream>
#include <stdio.h>
#include <fstream>
#include <math.h>
#include <string.h>
#include <stdlib.h>
double *reconstruct_laplacian_pyramid(double **pyr, double *filter, int **dimsList, int nlev, int *subwindow, int *dimsout);
#endif