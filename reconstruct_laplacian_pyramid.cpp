#include "upsample.h"
#include "child_window.h"
#include "reconstruct_laplacian_pyramid.h"

double *reconstruct_laplacian_pyramid(double **pyr, double *filter, int **dimsList, int nlev, int *subwindow, int *dimsout){
	int r = dimsList[0][1];//rows
	int c = dimsList[0][0];//cols
	
	int **subwindow_all = (int **)calloc(nlev, sizeof(int*));//size of nlev rows by 4 cols
	for(int i =0; i<nlev; i++)
	subwindow_all[i] = (int *)calloc(4, sizeof(int));
	
	if(subwindow == NULL){
	subwindow_all[0][0] = 1;
	subwindow_all[0][1] = r;
	subwindow_all[0][2] = 1;
	subwindow_all[0][3] = c;
	}
	else{
	//wont come here, because subwindow is passed in as NULL
	} 

	for(int lev = 1; lev<nlev; lev++){
	child_window(subwindow_all[lev-1], subwindow_all[lev]);
	}
	//check subwindow_all
	for(int i =0; i<nlev; i++){
		cout<<"subwindow_all: ["<<i<<"]:";
		for(int j= 0; j<4; j++)
		cout<<" ["<<subwindow_all[i][j]<<"]";
	cout<<endl;
	}

	//% start with low pass residual
	int kdimsout[3] = {0, 0, 0};
	double *R = pyr[nlev-1];
	for(int lev = nlev-2; lev>=0; lev--){
	// % upsample, and add to current level
	cout<<"R_dims: ["<<dimsList[lev+1][0]<<" "<<dimsList[lev+1][1]<<" "<<dimsList[lev+1][2]<<endl;	
	double *k = upsample(R, filter, dimsList[lev+1], kdimsout, subwindow_all[lev]);		
	cout<<"K_dimsout: ["<<kdimsout[0]<<" "<<kdimsout[1]<<" "<<kdimsout[2]<<endl;
		if(lev!=nlev-2)
		free(R);
	R = (double *)calloc(kdimsout[0]*kdimsout[1]*kdimsout[2], sizeof(double));
		for(int i =0; i<kdimsout[0]*kdimsout[1]*kdimsout[2]; i++)
		R[i]= pyr[lev][i] + k[i];	
	}
	dimsout[0] = kdimsout[0];
	dimsout[1] = kdimsout[1];
	dimsout[2] = kdimsout[2];
	return R;
}
