#include "upsample.h"
/***this function produces a matrix that is to be subtracted from J
****there is no dims change in this function either!
***/

double *upsample(double *I, double *filter, int *dims, int *dimsout, int *subwindow){
		
	//cout<<"int UPsample top"<<endl;
	//create filter matrix
	CImg<double> filtertemp (filter, 5, 5, 1, 1, false);
	//initialize the rest
	int r = subwindow[1] - subwindow[0] + 1;
	int c = subwindow[3] - subwindow[2] + 1;
	int k = dims[2]; 	//=3
	dimsout[0] = c;
	dimsout[1] = r;
	dimsout[2] = k;
	//cout<<"in upsample rck: "<<r<<" "<<c<<" "<<k<<" dimsout: ["<<dimsout[1]<<" "<<dimsout[0]<<" "<<dimsout[2]<<"]"<<endl;
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

	//cout<<"int UPsample 1"<<endl;
	//interpolate, convolve with 2D separable filter
	CImg<double> R(c, r, 1, k, 0);	        //initialize to 0 for R

	//cout<<"int UPsample 2"<<endl;
	int indexI = 0;
	for(int s= 0; s<k; s++){
		for(int i= reven; i<r; i=i+2){
			for(int j = ceven; j<c; j= j+2){
			//cout<<"j, i, s, index: "<<j<<" "<<i<<" "<<s<<" "<<indexI<<endl;
			R(j, i, 0, s) = I[indexI];
			indexI++;
			}
		} 
	}
	/*
	if(r==5&&c==5)
	R.save("./data/R_upsample.txt", -1);
	if(r==5&&c==6)
	R.save("./data/R_upsample1.txt", -1);
	if(r==5&&c==8)
	R.save("./data/R_upsample2.txt", -1);
	*/
/*	
	cout<<"check values1: "<<tempR[0]<<" "<<tempR[1]<<" "<<tempR[2]<<" "<<tempR[3]<<" "<<tempR[4]<<endl;
	cout<<"check values2: "<<tempR[dimsout[0]*dimsout[1]-1]<<" "<<tempR[dimsout[0]*dimsout[1]-2]<<" "<<tempR[dimsout[0]*dimsout[1]-3]<<endl;
	cout<<"check values3: "<<tempR[dimsout[0]*dimsout[1]*dimsout[2]-1]<<" "<<tempR[dimsout[2]*dimsout[0]*dimsout[1]-2]<<" "<<tempR[dimsout[2]*dimsout[0]*dimsout[1]-3]<<endl;
	*/
        R.correlate(filtertemp, 0, false);

	
	CImg<double> Z(c, r, 1, k, 0);
	for(int s= 0; s<k; s++){
		for(int i= reven; i<r; i=i+2){
			for(int j = ceven; j<c; j= j+2){
			Z(j, i, 0, s) = 1;
			}		
		}
	}
	/*if(r==5&&c==5)
	Z.save("./data/Z_upsample.txt", -1);
	if(r==5&&c==6)
	Z.save("./data/Z_upsample1.txt", -1);
	if(r==5&&c==8)
	Z.save("./data/Z_upsample2.txt", -1);	
	*/
	Z.correlate(filtertemp, 0, false);

	
	double *tempR = R;
	double *tempZ = Z;
	double *finalR= (double *)calloc(r*c*k, sizeof(double));
	for(int i = 0; i< r*c*k; i++){
	finalR[i] = tempR[i]/tempZ[i];
	}

	//cout<<"in UPsample bot"<<endl;
	return finalR;
}
