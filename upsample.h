#ifndef UPSAMPLE_H
#define UPSAMPLE_H
using namespace std;
#include <iostream>
#include <stdio.h>
#include <fstream>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include "CImg.h"
using namespace cimg_library;
//keep the dims as the input and dimsout is kept just in case it is needed in reconstruct.
double *upsample(double *I, double *filter, int *dims, int *dimsout, int *subwindow);//does not return subwindow_child
#endif
