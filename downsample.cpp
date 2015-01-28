#include "downsample.h"
#include "child_window.h"
double *downsample(double *I, double *filter, int *dims, int *dimsout, int *subwindow, int *subwindow_child){
	if(subwindow_child!=NULL){
	child_window(subwindow, subwindow_child);
	}
	//cout<<"int downsample top"<<endl;
	//low pass, convolve with 2D separable filter: first creat nonshare cimg image instance from I
	CImg<double> R(I ,dims[0], dims[1], 1, dims[2], false);	//for R
	//cout<<"int downsample 1"<<endl;
        //R = imfilter(I,filter);
	CImg<double> filtertemp (filter, 5, 5, 1, 1, false);
	//cout<<"int downsample 2"<<endl;
        R.correlate(filtertemp, 0, false);

	CImg<double> Z(dims[0], dims[1], 1, dims[2], 1);
	Z.correlate(filtertemp, 0, false);
	double *tempR = R;
	double *tempZ = Z;
	for(int i = 0; i< dims[0]*dims[1]*dims[2]; i++){
	tempR[i] = tempR[i]/tempZ[i];
	}
	//cout<<"int downsample 3"<<endl;
	
	//decimate
	int reven;
	int ceven;
	if(subwindow[0]%2==0)
	reven = 1;
	else 
	reven =  0;

	if(subwindow[2]%2==0)
	ceven = 1;
	else 
	ceven = 0;
	//cout<<"sub[0] sub[2]: ["<<subwindow[0]<<" "<<subwindow[2]<<"],"<<"rceven: ["<<reven<<" "<<ceven<<endl;
	//dims = w, h, chan
	int r = dims[1];//height max
	int c = dims[0];//width max
	dims[0] = round(((double)dims[0]-ceven)/2.0);
	dims[1] = round(((double)dims[1]-reven)/2.0);
	dims[2] = dims[2];	
	
	dimsout[0] = dims[0];
	dimsout[1] = dims[1];
	dimsout[2] = dims[2];

	//cout<<"r c: "<<r<<" "<<c<<endl;
	//cout<<"Down sample dims: "<<dims[0]<<" "<<dims[1]<<" "<<dims[2]<<endl;

	double *finalR= (double *)calloc(dims[0]*dims[1]*dims[2], sizeof(double));
	int index = 0;
	for(int s= 0; s<dims[2]; s++){
		for(int i= reven; i<r; i=i+2){
			for(int j = ceven; j<c; j= j+2){
			//cout<<"J I S index: "<<j<<" "<<i<<" "<<s<<" "<<index<<" "<<R(j, i, 0, s)<<endl;
			finalR[index] = R(j, i, 0, s);
			index++;
			}
		
		}
	}

	//cout<<"int downsample bot"<<endl;
	//cout<<"dims: "<<dims[0]<<" "<<dims[1]<<" "<<dims[2]<<" "<<dims[0]*dims[1]*dims[2]<<endl;
	return finalR;
}
