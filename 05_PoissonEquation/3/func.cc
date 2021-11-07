#include <iostream>
#include <cmath>
#include <vector>

using namespace std;

extern double p0;             
extern int lc;
extern int px;
extern int py;
extern int lx;
extern int ly;
extern int nx;
extern int ny;
extern double hx;
extern double hy;
extern bool debug;

// return l from i,j
double getl( int i, int j ){
	//if( debug ) std::cout<<(j-1)*nx+i-1<<std::endl;
	return (j-1)*nx+i-1;
}

// check if i,j is in non zero charge region
bool inLcP1( int i, int j){
	if(((i>=px/hx)&&(i<px/hx+0.5*lc/hx))&&((j>=py/hy)&&(j<py/hy+0.5*lc/hy))) return true;
	return false;
}
bool inLcN1( int i, int j){
	if(((i>=px/hx)&&(i<px/hx+0.5*lc/hx))&&((j>=py/hy+0.5*lc/hy)&&(j<py/hy+lc/hy))) return true;
	return false;
}
bool inLcP2( int i, int j){
	if(((i>=px/hx+0.5*lc/hx)&&(i<px/hx+lc/hx))&&((j>=py/hy+0.5*lc/hy)&&(j<py/hy+lc/hy))) return true;
	return false;
}
bool inLcN2( int i, int j){
	if(((i>=px/hx+0.5*lc/hx)&&(i<px/hx+lc/hx))&&((j>=py/hy)&&(j<py/hy+0.5*lc/hy))) return true;
	return false;
}

// return j from l
double getj( int l ){
	return l/nx+1;
}

// return i from l
double geti( int l ){
	return l-(getj(l)-1)*nx+1;
}

// return matrix element (l,m)
double getMlm( int l, int m ){
	// check if is an edge point
	if((geti(l)==1)||(geti(l)==nx)||(getj(l)==1)||(getj(l)==ny)){
		if(l==m) return 1;
		else return 0;
	}
	
	if(((geti(l)==geti(m)-1)&&(getj(l)==getj(m)))||((geti(l)==geti(m)+1)&&(getj(l)==getj(m)))) return pow(hx,-2);
	if(((getj(l)==getj(m)-1)&&(geti(l)==geti(m)))||((getj(l)==getj(m)+1)&&(geti(l)==geti(m)))) return pow(hy,-2);
	if(l==m) return -2.*(pow(hx,-2)+pow(hy,-2));
	return 0;
}

// compute difference between iterations
double computediff( vector<double> phi_0, vector<double> phi_1 ){
	double sumdiff = 0;
	for(int l = 0; l< nx*ny; ++l) sumdiff+= abs(phi_0.at(l)-phi_1.at(l))/(nx*ny);
	return sumdiff;
}