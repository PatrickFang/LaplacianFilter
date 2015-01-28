#ifndef DOWNSAMPLE_H
#define DOWNSAMPLE_H
using namespace std;
#include <iostream>
#include <stdio.h>
#include <fstream>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include "CImg.h"
using namespace cimg_library;
double* downsample(double *I, double *filter, int *dims, int *dimsout, int *subwindow, int *subwindow_child);//return subwindow_child too
#endif
