#include "child_window.h"
void child_window(int *parent, int *child){
	int N = 1;
	
	//cout<<"parent printed: "<<parent[0]<<" "<<parent[1]<<" "<<parent[2]<<" "<<parent[3]<<endl;
	for (int i =0; i<N; i++){
	child[0] = ceil((double)(parent[0]+1)/2.0);
	child[2] = ceil((double)(parent[2]+1)/2.0);
	child[1] = floor((double)(parent[1]+1)/2.0);
	child[3] = floor((double)(parent[3]+1)/2.0);
	}	
} 
