#include "numlevels.h"
#include "downsample.h"
#include "upsample.h"
#include "laplacian_pyramid.h" 

double **laplacian_pyramid(double *I, double *filter, int *dims, int nlev, int *subwindow, int **dimsout){
	//I comes as zeros(mxnx3)
	//cout<<"lp___pyramid top"<<endl;
	int r = dims[1];
	int c = dims[0];
	if(nlev == -1){
	int im_sz[2]  = {r, c};	
	int nlev = numlevels(im_sz);
	}
	if(subwindow ==NULL)
	int subwindow[4] = {1, r, 1, c};


	//declare L container	
	double **L = (double **)calloc(nlev, sizeof(double *));
	int dimsList_dummy[3] = {0, 0, 0};	//this variable is no use for now
						//assuming that the laplacian pyramid will
						//have exactly the same dimensions as the G_pyramid
	
	
	double *J = I;
	//recursively call downsample, upsample and let J = I; build pyramid
	for(int i = 0; i<nlev-1; i++){ 
	//cout<<"lp___pyramid middle"<<endl;	
	//cout<<"subwindow printed: "<<subwindow[0]<<" "<<subwindow[1]<<" "<<subwindow[2]<<" "<<subwindow[3]<<endl;
	//apply low pass filter, and downsample
	dimsout[i] = (int *)calloc(3, sizeof(int));
	dimsout[i][0] = dims[0];
	dimsout[i][1] = dims[1];
	dimsout[i][2] = dims[2];
	int subwindow_child[4] = {0, 0, 0 ,0};
	double *tempI = I;
	//cout<<"lp__pyramid dims before downsample: "<<dims[0]<<" "<<dims[1]<<" "<<dims[2]<<"______________"<<endl;
	//cout<<"lp__pyramid dims_out before downsample!: "<<dimsout[i][0]<<" "<<dimsout[i][1]<<" "<<dimsout[i][2]<<"_______________"<<endl;
	
	I = downsample(J, filter, dims, dimsList_dummy, subwindow, subwindow_child);
	//cout<<"lp__pyramid dims after downsample: "<<dims[0]<<" "<<dims[1]<<" "<<dims[2]<<"_______________"<<endl;
	//cout<<"lp__pyramid dims_out after downsample!: "<<dimsout[i][0]<<" "<<dimsout[i][1]<<" "<<dimsout[i][2]<<"_______________"<<endl;
	
	//in each level, store difference between image and upsampled low pass version	
	double *upsample_I = upsample(I, filter, dims, dimsList_dummy, subwindow);
	/*cout<<"lp__pyramid dims after UPsample!: "<<dims[0]<<" "<<dims[1]<<" "<<dims[2]<<"_______________"<<endl;
	cout<<"lp__pyramid dims_out after UPsample!: "<<dimsout[i][0]<<" "<<dimsout[i][1]<<" "<<dimsout[i][2]<<"_______________"<<endl;
	cout<<"lp__pyramid dims_dummy after UPsample!: "<<dimsList_dummy[0]<<" "<<dimsList_dummy[1]<<" "<<dimsList_dummy[2]<<"_______________"<<endl;
	
	cout<<"check values1: "<<upsample_I[0]<<" "<<upsample_I[1]<<" "<<upsample_I[2]<<" "<<upsample_I[3]<<" "<<upsample_I[4]<<endl;
	cout<<"check values2: "<<upsample_I[dimsList_dummy[0]*dimsList_dummy[1]-1]<<" "<<upsample_I[dimsList_dummy[0]*dimsList_dummy[1]-2]<<" "<<upsample_I[dimsList_dummy[0]*dimsList_dummy[1]-3]<<endl;
	cout<<"check values3: "<<upsample_I[dimsList_dummy[0]*dimsList_dummy[1]*dimsList_dummy[2]-1]<<" "<<upsample_I[dimsList_dummy[2]*dimsList_dummy[0]*dimsList_dummy[1]-2]<<" "<<upsample_I[dimsList_dummy[2]*dimsList_dummy[0]*dimsList_dummy[1]-3]<<endl;
	*/
	

	L[i] = (double *)calloc(dimsout[i][0]*dimsout[i][1]*dimsout[i][2], sizeof(double));
		for(int j =0; j<dimsout[i][0]*dimsout[i][1]*dimsout[i][2]; j++)
		L[i][j] = J[j] - upsample_I[j];
	free(upsample_I);
	free(tempI);
	
	//continue with low pass image
	J = I;
		//copy subwindow_child into subwindow
		for (int k = 0; k<4; k++)
		subwindow[k] = subwindow_child[k];
		
	}

	L[nlev-1] = J;
	dimsout[nlev-1] = (int *)calloc(3, sizeof(int));
	dimsout[nlev-1][0] = dims[0];
	dimsout[nlev-1][1] = dims[1];
	dimsout[nlev-1][2] = dims[2];
	//cout<<"lp___pyramid bot"<<endl;
	return L;
}
