using namespace std;
#include <iostream>
#include <stdio.h>
#include <fstream>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include "numlevels.h"
int numlevels(int *im_sz){
	int nlevels= 1;
	int min_d = min(im_sz[0], im_sz[1]);
	while (min_d >1){
	nlevels++;
	min_d = floor((min_d+1)/2);
	}
	return nlevels;
}
