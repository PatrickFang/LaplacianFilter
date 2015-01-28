using namespace std;
#include <iostream>
#include <stdio.h>
#include <fstream>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include "CImg.h"
#include "numlevels.h"
#include "downsample.h"
#include "laplacian_pyramid.h"
#include "upsample.h"
#include "reconstruct_laplacian_pyramid.h"
using namespace cimg_library;
/****some helper functions****/
static inline double min(double x, double y) { return (x <= y ? x : y); }	//function overloading 
static inline double max(double x, double y) { return (x <= y ? y : x); }

static inline int min(int x, int y) { return (x <= y ? x : y); }
static inline int max(int x, int y) { return (x <= y ? y : x); }
#define COLOR_REMAPPING 1	//if 1==>rbg, else lum
#define DOMAIN_LINLOG 1		//if 1==>lin, else log
#define DISPLAY_WINDOW 1	//if 1==>enable display, else no displaying
#define RESIZE_ENABLE 1
/*

% convert RGB to grayscale intensity
function Y = luminance(I)
    switch size(I,3),
        case 1, Y = I;
        case 3, Y = (20*I(:,:,1) + 40*I(:,:,2) + I(:,:,3))/61;
    end
end



% grayscale remapping function
function inew = r_gray(i,g0,sigma_r,fd,fe)
    dnrm = abs(i-g0);
    dsgn = sign(i-g0);
    % detail and edge processing
    rd = g0 + dsgn*sigma_r.*fd(dnrm/sigma_r);
    re = g0 + dsgn.*(fe(dnrm - sigma_r) + sigma_r);
    % edge-detail separation based on sigma_r threshold
    isedge = dnrm > sigma_r;
    inew = ~isedge.*rd + isedge.*re;
end


% color remapping function
function inew = r_color(i,g0,sigma_r,fd,fe)
    g0 = repmat(g0,[size(i,1) size(i,2) 1]);
    dnrm = sqrt(sum((i-g0).^2,3));
    unit = (i-g0)./repmat(eps + dnrm,[1 1 3]);
    % detail and edge processing
    rd = g0 + unit.*repmat(sigma_r*fd(dnrm/sigma_r),[1 1 3]);
    re = g0 + unit.*repmat((fe(dnrm - sigma_r) + sigma_r),[1 1 3]);
    % edge-detail separation based on sigma_r threshold
    isedge = repmat(dnrm > sigma_r,[1 1 3]);
    inew = ~isedge.*rd + isedge.*re;
end
*/





/*
 edge remapping function
*/
double *fe(double *a, int *dims, double beta, double sigma_r){
	double *out = (double *)calloc(dims[0]*dims[1]*dims[2], sizeof(double));
	for(int i = 0; i<dims[0]*dims[1]*dims[2];i++){
	out[i] = (a[i] - sigma_r)*beta;
	}
	return out;
}
/*** 
smooth step edge between (xmin,0) and (xmax,1)
***/
double *smooth_step(double xmin, double xmax, double *x, int *dims){
	double *y = (double *)calloc(dims[0]*dims[1]*dims[2], sizeof(double));
	for(int i = 0; i<dims[0]*dims[1]*dims[2]; i++){
	y[i] = (x[i] - xmin)/(xmax - xmin);
		if(y[i] > 1)	
		y[i] = 1;
		else if(y[i] <0)
		y[i] = 0;	
	y[i] = pow(y[i], 2)*pow(y[i]-2, 2);
	}	
	return y;
}

/**
input: alpha, level,sigma_r,d =matrix, and dims of d
help functions: smooth_step
//takes 2 arguements: a matrix dnrm/sigma_r and alpha
**/
double *fd(double *d, int *dims, double alpha, double noise_level, double sigma_r){
	//the dims[2] = 1 when it is called from r_color or it might be 3
	double *out = (double *)calloc(dims[0]*dims[1]*dims[2], sizeof(double));
	for (int i =0; i<dims[0]*dims[1]*dims[2]; i++){
	out[i] = pow(d[i], alpha);
	}
	if (alpha<1){//call smooth step func. to produce tau
	double *dxs_r = (double *)calloc(dims[0]*dims[1]*dims[2], sizeof(double));
		for(int i = 0; i<dims[0]*dims[1]*dims[2]; i++)
		dxs_r[i] = d[i]*sigma_r;

	double *tau = smooth_step(noise_level, 2*noise_level, dxs_r, dims);
		for(int i = 0; i<dims[0]*dims[1]*dims[2]; i++)
		out[i] = tau[i]*out[i] + (1-tau[i])*d[i];
	free(dxs_r);
	free(tau);
	}	
	return out;
}



/**
input: sigma_r, g0 and isub
help functions: fd and fe
this function remaps g0 into isub's dimensions
and manipulate them to produce Iremap of Isub_dims==>dims
**/
#if COLOR_REMAPPING
double *r_function(double sigma_r, double *g0temp, double *Isub, int *dims, double alpha, double noise_level, double beta){
	//cout<<"Isub dims: ["<<dims[0]<<" "<<dims[1]<<" "<<dims[2]<<"]"<<endl;
	//first remap g0temp to g0 of Isub
	double *g0 = (double *)calloc(dims[0]*dims[1]*dims[2], sizeof(double));
	double **dnrmtemp = (double **)calloc(dims[2], sizeof(double *));//used to temporaly record values of 
								   //dnrm sum and sqrt
	double *dnrm = (double *)calloc(dims[0]*dims[1], sizeof(double));
	for(int s=0; s<dims[2]; s++){
	dnrmtemp[s] = (double *)calloc(dims[0]*dims[1], sizeof(double)); 
		for(int i =0; i<dims[0]*dims[1]; i++){
		g0[i + s*dims[0]*dims[1]] = g0temp[s];
		dnrmtemp[s][i] = pow((Isub[i + s*dims[0]*dims[1]] - g0[i + s*dims[0]*dims[1]]), 2);
		}
	}
	
	//sum up along the three dims and sqrt dnrmtemp to become dnrm
	//all produce unit = (i-g0)./repmat(eps + dnrm,[1 1 3]);
	double *unit = (double *)calloc(dims[0]*dims[1]*dims[2], sizeof(double));
	for(int i= 0; i<dims[0]*dims[1]; i++){
		double t = 0;
		for(int s=0; s<dims[2]; s++){
			//dnrm[i] = sqrt(dnrmtemp[0][i] + dnrmtemp[1][i] + dnrmtemp[2][i]);
			t = t+dnrmtemp[s][i];
			}
		dnrm[i] = sqrt(t);
		
		if(dnrm[i] == 0){
			for(int s=0; s<dims[2]; s++){
			//unit[i] = 0;
			//unit[i + dims[0]*dims[1]] = 0;
			unit[i + dims[0]*dims[1]*s] = 0;
			}
		}
		else{
			for(int s=0; s<dims[2]; s++){
			//unit[i] = (Isub[i] - g0[i])/dnrm[i];
			//unit[i + dims[0]*dims[1]] = (Isub[i + dims[0]*dims[1]] - g0[i + dims[0]*dims[1]])/dnrm[i]; 
			unit[i + dims[0]*dims[1]*s] = (Isub[i + s*dims[0]*dims[1]] - g0[i + s*dims[0]*dims[1]])/dnrm[i];
			}
		}
		
	}

	// detail and edge processing
	int fddims[3] = {dims[0], dims[1], 1};
	double *dnrmscaled = (double *)calloc(fddims[0]*fddims[1]*fddims[2], sizeof(double));
		for(int i = 0; i<fddims[0]*fddims[1]*fddims[2]; i++)
		dnrmscaled[i] = dnrm[i]/sigma_r;
	double *fd_result = fd(dnrmscaled, fddims, alpha, noise_level, sigma_r);
	double *fe_result = fe(dnrm, fddims, beta, sigma_r);
	double *rd  = (double *)calloc(dims[0]*dims[1]*dims[2], sizeof(double));
	double *re  = (double *)calloc(dims[0]*dims[1]*dims[2], sizeof(double));
	int fd_index = 0;
	for(int i = 0; i<dims[0]*dims[1]*dims[2]; i++){
		if(i%(dims[0]*dims[1])==0)
		fd_index = 0;
	rd[i] = g0[i] + unit[i]*sigma_r*fd_result[fd_index];
	re[i] = g0[i] + unit[i]*(fe_result[fd_index] + sigma_r);
	fd_index++;
	}
	
	// edge-detail separation based on sigma_r threshold
	double *isedge  = (double *)calloc(dims[0]*dims[1]*dims[2], sizeof(double));
	double *isedge_inverse = (double *)calloc(dims[0]*dims[1]*dims[2], sizeof(double));
	double *inew  = (double *)calloc(dims[0]*dims[1]*dims[2], sizeof(double));
	//isedge is a matrix of values of 1s and 0s
	int index_dnrm = 0;
	//this loop produce isedge as well as the inew  => Iremap
    	for(int i = 0; i<dims[0]*dims[1]*dims[2]; i++){
		if(i%(dims[0]*dims[1])==0)
		index_dnrm = 0;
		if(dnrm[index_dnrm]>sigma_r){
		isedge[i] = 1;
		isedge_inverse[i] = 0;
		}
		else{ 
		isedge[i] = 0;
		isedge_inverse[i]= 1;
		}
	index_dnrm++;
	inew[i] = isedge[i]*re[i] + (isedge_inverse[i]*rd[i]);
	}
	/*			
	//chck g0 now
	for(int i = 0; i<dims[0]*dims[1]*dims[2]; i++){
	cout<<g0[i]<<" ";
	}
	cout<<endl;
	*/
	//free everything
	free(g0);
	free(dnrm);
	free(dnrmscaled);
	free(fe_result);
	free(fd_result);
	for(int i= 0; i<3; i++)
	free(dnrmtemp[i]);
	free(dnrmtemp);
	free(unit);	
	free(rd);
	free(re);
	free(isedge);	
	free(isedge_inverse);
	return inew;
}
#else //lum mode




#endif


int main(int argc, char* argv[]){
	if(argc<2){
	cout<<"error: input the image name"<<endl;
	return 0;
	}
	//this part will later become arguements of main
	double sigma_r; 
	double alpha; 
	double beta; 
	char *stopstring= NULL;
	
	if(argc ==2){
	cout<<"Warning: missing arguments sigma_r, alpha, beta."<<endl;
	cout<<"Using default values:\nsigma_r = 0.4\nalpha = 0.25\nbeta = 1"<<endl;
	sigma_r = 0.4;
	alpha = 0.25;
	beta = 1;
	}
	
	else if(argc ==3){	
	sigma_r = strtod(argv[2], &stopstring);
		if(stopstring[0]!='\0'){
		cout<<"error: parameters are double"<<endl;
		return -1;
		}
	cout<<"Warning: missing arguments alpha, beta."<<endl;
	cout<<"Using default values:\nalpha = 0.25\nbeta = 1"<<endl;
	alpha = 0.25;
	beta = 1;	
	}

	else if(argc ==4){
	sigma_r = strtod(argv[2], &stopstring);
		if(stopstring[0]!='\0'){
		cout<<"error: parameters are double"<<endl;
		return -1;
		}
	stopstring = NULL;
	alpha = strtod(argv[3], &stopstring);
		if(stopstring[0]!='\0'){
		cout<<"error: parameters are double"<<endl;
		return -1;	
		}
	cout<<"Warning: missing arguments beta."<<endl;
	cout<<"Using default values:\nbeta = 1"<<endl;
	beta = 1;
	}

	else if(argc==5){
	sigma_r = strtod(argv[2], &stopstring);
		if(stopstring[0]!='\0'){
		cout<<"error: parameters are double"<<endl;
		return -1;
		}
	stopstring = NULL;
	alpha = strtod(argv[3], &stopstring);
		if(stopstring[0]!='\0'){
		cout<<"error: parameters are double"<<endl;
		return -1;	
		}
	stopstring = NULL;
 	beta = strtod(argv[4], &stopstring);
		if(stopstring[0]!='\0'){
		cout<<"error: parameters are double"<<endl;
		return -1;	
		}
	}

	cout<<"parameters: "<<sigma_r<<" "<<alpha<<" "<<beta<<endl;
	char *imagename = argv[1];
	char imagepath[50];
	strcpy(imagepath, "");
	strcat(imagepath, imagename);
	cout<<imagepath<<endl;
	//read in the image	
	CImg<double>image(imagepath);
	
	//image.save("./data/cppimage.txt", -1);

	//some specfications of the image
	int width = image.width();		//cols
	int height = image.height();		//rows
	int depth = image.depth(); 		//depth
	int spec = image.spectrum();
	if(spec!=3){
	cout<<"Error: The input image has to be RBG and contains three color channels"<<endl;
	return -1;
	}
	int totalnumpixels = image.size();
	cout<<height<<" "<<width<<" "<<depth<<" "<<spec<<" "<<totalnumpixels<<endl;
	double *temp = image;
	for(int i = 0; i< totalnumpixels; i++){
	temp[i] = temp[i]/255;
	}
	
	//resize by Bicubic interpolation: optional here
	#if RESIZE_ENABLE
	image.resize(ceil((double)width/4.0), ceil((double)height/4.0), 1, 3, 5, 0);
	#endif 
	int h = image.height();
	int w = image.width();
	int chan = image.spectrum();
	//make sure the pixels have values less than 1 but greater than 0
	temp =image;
	for(int i = 0; i<h*w*spec; i++){
		if(temp[i]<0)
		temp[i] = 0;
		if(temp[i]>1)
		temp[i] = 1;
	}
	cout<<"new dims:" <<h<<" "<<w<<" "<<spec<<endl;
	//image.save("./data/scaledimage255.txt", -1);

	char *colorRemapping = "rgb"; //could be "lum"
	
	//NOW: call R = lapfilter(I,sigma_r,alpha,beta,colorRemapping,domain);

	#if DOMAIN_LINLOG == 0	//log
		sigma_r = log(sigma_r);
	#endif

	//detailed remapping function
	double noise_level = 0.01;
	
	/*
	if(strcmp(colorRemapping, "rgb") == 0){
	//define function r as
	//r = @(i,g0)(r_color(i,g0,sigma_r,@fd,@fe));	
	}
	else if(strcmp(colorRemapping, "lum") == 0){
	//define function r as
	//r = @(i,g0)(r_gray(i,g0,sigma_r,@fd,@fe));
	
	//and manipulate I
	}
	else 
	cout<<"error: invalid color remapping"<<endl;
	*/
	



	//if(alpha==1&&beta==1)
	//R = I; do nothing
	
	//NOW call lapfilter_core	
	//first wanna create G pyramid
	//declare filter 
	double filter[25]={0.0025,
			0.0125,
			0.02,
			0.0125,
			0.0025,
			0.0125,
			0.0625,
			0.1,
			0.0625,
			0.0125,
			0.02,
			0.1,
			0.16,
			0.1,
			0.02,
			0.0125,
			0.0625,
			0.1,
			0.0625,
			0.0125,
			0.0025,
			0.0125,
			0.02,
			0.0125,
			0.0025
			};	


	//get the image txt file from matlab
	/*******IMPORTANT FOR TESTING ****/
	/*
	CImg <double>matlabI("./data/IbeforeIsub.txt");	//flower	
	//CImg <double>matlabI("./data/IbeforeIsub_beach.txt");	
	double *matlab = matlabI;
	CImg <double>matremap(w, h, 1, spec);
	int index_re = 0;
	int x0 =0, x1=0, x2=0, y0=0, y1=0, y2=0;
		while(y2<h){
			
			int multipleof200 = index_re/w;
			if(multipleof200%3 == 0){
			matremap(x0, y0, 0, 0) = matlab[index_re];
				x0++;
				if(x0 ==w){
				x0 =0;
				y0++;
				}
			}
			else if(multipleof200%3 == 1){
			matremap(x1, y1, 0, 1) = matlab[index_re]; 
			x1++;
				if(x1==w){
				x1 =0;
				y1 ++;
				}
			}
			else if(multipleof200%3 == 2){
			matremap(x2, y2, 0, 2) = matlab[index_re]; 
				x2++;
				if(x2 ==w){
				x2 =0;
				y2 ++;
				}

			}
			index_re++;
		}
	matlab = matremap;
		

	#if DOMAIN_LINLOG == 0 	//log mode: gotta do to_main function
			//if is 'lin' mode, do nothing!
        //to_domain   = @(I) log(I + eps);
	for(int i =0; i<h*w*spec; i++){
		if(matlab[i]==0)
		matlab[i] = -36.044;	//in matlab log(0+eps) =-36.044
		else 
		matlab[i]  = log(matlab[i]);
	}
	matremap.save("./data/matlabI_log_rbg.txt", -1);	
	#endif 
	*/
	/*
	cout<<endl;
	cout<<"check matlab1: "<<matlab[0]<<" "<<matlab[1]<<" "<<matlab[2]<<endl;
	cout<<"check matlabmiddle: "<<matlab[199]<<" "<<matlab[200]<<" "<<matlab[201]<<endl;
	cout<<"check matlabmiddle2: "<<matlab[599]<<" "<<matlab[600]<<" "<<matlab[601]<<endl;
	cout<<"check matlab2: "<<matlab[w*h*spec-1]<<" "<<matlab[w*h*spec-2]<<" "<<matlab[w*h*spec-3]<<endl;
	*/
	/*******IMPORTANT FOR TESTING ENDS****/
	//cout<<"//CALL GAUSSIAN_PYRAMID"<<endl;
	//CALL GAUSSIAN_PYRAMID
	int r = h;
	int c = w;
	int im_sz[2]  = {r, c};		
	int nlev = numlevels(im_sz);
	cout<<nlev<<endl;
	int subwindow[4] = {1, r, 1, c};
	//declare G container
	double **G = (double **)calloc(nlev, sizeof(double *));
	int **dimsList = (int **)calloc(nlev, sizeof(int *));
	//first cell is just the image, here maybe want to use a non-share!
	G[0] = image;
	//G[0] = matremap;	//!!!!!!!!!!!!!!for testing purpose here i used matlab image not C image!\
				//NOTE the result is much better!	
	//cout<<"//CALL GAUSSIAN_PYRAMID2"<<endl;
	dimsList[0] = (int *)calloc(3, sizeof(int));
	dimsList[0][0] = w;
	dimsList[0][1] = h;
	dimsList[0][2] = spec;
	int dims[3] = {w, h, spec};
	//call downsample
	for(int i= 1; i<nlev; i++){
	//recursively call downsample with previous level of G
	dimsList[i] = (int *)calloc(3, sizeof(int));
	G[i] = downsample(G[i-1], filter, dims, dimsList[i], subwindow, NULL);
	}
	cout<<"subwindow G: ["<<subwindow[0]<<" "<<subwindow[1]<<" "<<subwindow[2]<<" "<<subwindow[3]<<endl;

	//check
	/*
	for (int i = 0; i<nlev; i++)
	{
	cout<<"dimsList: "<<dimsList[i][0]<<" "<<dimsList[i][1]<<" "<<dimsList[i][2]<<endl;
	cout<<"G: "<<G[i][0]<<" "<<G[i][1]<<" "<<G[i][2]<<endl;
	cout<<"G: "<<G[i][dimsList[i][0]*dimsList[i][1]*dimsList[i][2]-1]<<" "<<G[i][dimsList[i][0]*dimsList[i][1]*dimsList[i][2]-2]<<" "<<G[i][dimsList[i][0]*dimsList[i][1]*dimsList[i][2]-3]<<endl;
	}
	*/
	
	
	//CALL LAPLACIAN_PYRAMID: produce pyramid of the same dims as G
	double *Izeros  = (double *)calloc(h*w*spec, sizeof(double));
		for(int i = 0; i<h*w*spec; i++)
		Izeros[i] = 0;
	dims[0] = w;
	dims[1] = h;
	dims[2] = spec;
	int **L_dims = (int **)calloc(nlev, sizeof(int *));
	double **L = laplacian_pyramid(Izeros, filter, dims, nlev, subwindow, L_dims);
	
	cout<<"subwindow L: ["<<subwindow[0]<<" "<<subwindow[1]<<" "<<subwindow[2]<<" "<<subwindow[3]<<endl;

	//check
	/*
	for (int i = 0; i<nlev; i++)
	{
	cout<<"L_dims: "<<L_dims[i][0]<<" "<<L_dims[i][1]<<" "<<L_dims[i][2]<<endl;
	cout<<"L: "<<L[i][0]<<" "<<L[i][1]<<" "<<L[i][2]<<" "<<L[i][3]<<" "<<L[i][4]<<" "<<L[i][5]<<endl;
	}
	*/


	//NOW the big loop in core
	for(int lev0 = 0; lev0<nlev; lev0++){
	int hw = 3*pow(2, lev0+1)-2;
	//cout<<"hw: "<<hw<<endl;
	cout<<"max dims: "<<dimsList[lev0][1]<<" "<<dimsList[lev0][0]<<endl;
	int flag = 0;
	//ignore the fprintf statements here
		for(int y0 = 0; y0<dimsList[lev0][1]; y0++){
			for(int x0 = 0; x0<dimsList[lev0][0]; x0++){
			
			// % coords in full-res image corresponding to (lev0,y0,x0)
			int yf = y0*pow(2, lev0);
			int xf = x0*pow(2, lev0);
			// % subwindow in full-res image needed to evaluate (lev0,y0,x0) in result
			int yrng[2] = {max(0, yf-hw), min(dimsList[0][1]-1, yf+hw)}; 
			int xrng[2] = {max(0, xf-hw), min(dimsList[0][0]-1, xf+hw)};

			
			//if(lev0 == 3&& y0==0){//extra starts
			//cout<<"new x0: "<<x0<<" max x0: "<<dimsList[lev0][0]<<endl;
			//cout<<"yf xf: ["<<yf<<" "<<xf<<"] "<<endl;; 	
			//cout<<"yrng: ["<<yrng[0]<<" "<<yrng[1]<<"] "; 	
			//cout<<"xrng: ["<<xrng[0]<<" "<<xrng[1]<<"] "<<endl; 	
			//flag= 1;
					
			
			//copy the sub set of I of size yrng by xrng
			int Isub_dims[3] = {(xrng[1]-xrng[0]+1), (yrng[1]-yrng[0]+1), spec}; //w, h, chan
			int Isubdims[3] = {(xrng[1]-xrng[0]+1), (yrng[1]-yrng[0]+1), spec};
			double *Isub = (double *)calloc((yrng[1]-yrng[0]+1)*(xrng[1]-xrng[0]+1)*spec, sizeof(double));
			int index=0;
			for(int s= 0; s<spec; s++){
				for(int imy= yrng[0]; imy<=yrng[1]; imy++){
					for(int imx = xrng[0]; imx<=xrng[1]; imx++){
					Isub[index] = image(imx, imy, 0, s);
					//Isub[index] = matremap(imx, imy, 0, s); //with matlab image
					index++;
					}
				}
			}
			// % use the corresponding Gaussian pyramid coefficient to remap
           		// % the full-res subwindow
			int offset = y0*dimsList[lev0][0]+x0;//offset of x, y coordinates.
			double g0temp[3];
			g0temp[0] = G[lev0][offset];
			if(spec==3){
			g0temp[1] = G[lev0][offset+dimsList[lev0][0]*dimsList[lev0][1]];
			g0temp[2] = G[lev0][offset+dimsList[lev0][0]*dimsList[lev0][1]*2];
			}
	/*
			cout<<"x0 y0: ["<<x0<<" "<<y0<<"] offset: "<<offset<<" "<<offset+dimsList[lev0][0]*dimsList[lev0][1]<<" "<<offset+dimsList[lev0][0]*dimsList[lev0][1]*2<<endl;
			*/
			//call function r_function:  r could r_color or r_gray.. now call it r_function
			double *Iremap  = r_function(sigma_r, g0temp, Isub, Isub_dims, alpha, noise_level, beta);	

			// compute Laplacian pyramid for remapped subwindow
			int subwindow_Lremap[4] = {yrng[0]+1, yrng[1]+1, xrng[0]+1, xrng[1]+1};
			int **Lremap_dims = (int **)calloc(lev0+2, sizeof(int *));
            		double **Lremap = laplacian_pyramid(Iremap, filter, Isub_dims, lev0+2, subwindow_Lremap, Lremap_dims);
            	

				//check Lremap
				/*
				for (int i = 0; i<lev0+2; i++)
				{
				cout<<"Lremap_dims["<<i<<"]: "<<Lremap_dims[i][0]<<" "<<Lremap_dims[i][1]<<" "<<Lremap_dims[i][2]<<endl;
				cout<<"Lremap: "<<Lremap[i][0]<<" "<<Lremap[i][1]<<" "<<Lremap[i][2]<<" "<<Lremap[i][3]<<" "<<Lremap[i][4]<<" "<<Lremap[i][5]<<endl;
				cout<<"Lremap2: "<<Lremap[i][Lremap_dims[i][0]*Lremap_dims[i][1]-1]<<" "<<Lremap[i][Lremap_dims[i][0]*Lremap_dims[i][1]-2]<<" "<<Lremap[i][Lremap_dims[i][0]*Lremap_dims[i][1]-3]<<endl;
				cout<<"Lremap3: "<<Lremap[i][Lremap_dims[i][0]*Lremap_dims[i][1]*3-1]<<" "<<Lremap[i][3*Lremap_dims[i][0]*Lremap_dims[i][1]-2]<<" "<<Lremap[i][3*Lremap_dims[i][0]*Lremap_dims[i][1]-3]<<endl;
				}
				*/			
			// % bookkeeping to compute index of (lev0,y0,x0) within the
            		// % subwindow, at full-res and at current pyramid level
			int yfc = yf - yrng[0];
           		int xfc = xf - xrng[0];
            		int yfclev0 = floor(yfc/pow(2,lev0));
            		int xfclev0 = floor(xfc/pow(2,lev0));
			//cout<<"yfc xfc: ["<<yfc<<" "<<xfc<<"] yxlev0: ["<<yfclev0<<" "<<xfclev0<<"]"<<endl;
			
			// % set coefficient in result based on the corresponding
           		// % coefficient in the remapped pyramid
			int off_fc = yfclev0*Lremap_dims[lev0][0]+xfclev0;
			L[lev0][offset] = Lremap[lev0][off_fc];       
			if(spec==3){
			L[lev0][offset+L_dims[lev0][0]*L_dims[lev0][1]] = Lremap[lev0][off_fc + Lremap_dims[lev0][0]*Lremap_dims[lev0][1]];   
			L[lev0][offset+L_dims[lev0][0]*L_dims[lev0][1]*2] = Lremap[lev0][off_fc + Lremap_dims[lev0][0]*Lremap_dims[lev0][1]*2];
			}
			//cout<<"L[lev0][off+..]: "<<L[lev0][offset]<<" "<<L[lev0][offset+L_dims[lev0][0]*L_dims[lev0][1]]<<" "<<L[lev0][offset+L_dims[lev0][0]*L_dims[lev0][1]*2]<<endl;
			
			free(Isub);
			for(int i = 0; i<lev0+2; i++){
			free(Lremap[i]);		
			free(Lremap_dims[i]);
			}
			free(Lremap);
			free(Lremap_dims);	
			
			//} //extra ends
	
				
	
					
			}//x0 loop ends

			//if(flag){
			//sleep(3);
			//flag=0;
			//}

		}//y0 loop ends
	}//lev0 loop ends
	
	cout<<"Lend: "<<L[nlev-1][0]<<" "<<L[nlev-1][1]<<" "<<L[nlev-1][2]<<" "<<L[nlev-1][3]<<" "<<L[nlev-1][4]<<" "<<L[nlev-1][5]<<endl;
	cout<<"Gend: "<<G[nlev-1][0]<<" "<<G[nlev-1][1]<<" "<<G[nlev-1][2]<<" "<<G[nlev-1][3]<<" "<<G[nlev-1][4]<<" "<<G[nlev-1][5]<<endl;
	free(L[nlev-1]);	
	L[nlev-1] = G[nlev-1];  //% residual not affected
	//CALL reconstruct_l_p
	int R_dimsout[3] = {0, 0, 0};
	double *R = reconstruct_laplacian_pyramid(L, filter, dimsList, nlev, NULL, R_dimsout);

	//% postprocessing
	#if DOMAIN_LINLOG == 0	//only compiles if log mode
	//first do from_domain function:  @(R) exp(R) - eps;
	for(int i= 0; i<w*h*spec; i++){
	}
	CImg <double>lala1(R, R_dimsout[0], R_dimsout[1], 1, 3, false);
	lala1.save("./data/Rfinal_flower_log_rbg_withfromdomain.txt", -1);	
	//here the real postprocessing starts
	if(beta<=1){		
	}
	#endif 
	


	#if COLOR_REMAPPING == 0
	
	#endif



	//% clip out of bounds intensities
	for (int i= 0; i<h*w*spec; i++){
		if(R[i]<0)
		R[i] = 0;
		if(beta<1)
			if(R[i]>1)
			R[i] = 1;
	}


	#if DOMAIN_LINLOG == 0
	double gamma_val = 2.2;
	for (int i= 0; i < w*h*spec; i++)	
	R[i] = pow(R[i], 1/gamma_val);
	#endif 
	
	

	CImg <double>Rimage(R, R_dimsout[0], R_dimsout[1], 1, 3, false);
	//Rimage.save("./data/Rfinal_beachmatlab.txt", -1);	
	#if DISPLAY_WINDOW
	image.display();
	Rimage.display();
	#endif 
	//matlabI.display();
	//free dimsList and Gpyramid
	for(int i=0; i<nlev; i++){
	free(dimsList[i]);
		if(i!=nlev-1&&i!=0)
		free(G[i]);
	free(L[i]);
	free(L_dims[i]);
	}
	free(L_dims);
	free(L);
	free(G);
	free(dimsList);

	return 0;
}


