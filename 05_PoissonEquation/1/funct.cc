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



// check if i,j is in non zero charge region
bool inLc( int i, int j){
	if(((i>=px/hx)&&(i<px/hx+lc/hx))&&((j>=py/hy)&&(j<py/hy+lc/hy))) return true;
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

// return l from i,j
double getl( int i, int j ){
	//if( debug ) std::cout<<(j-1)*nx+i-1<<std::endl;
	double l=(j-1)*nx+i-1;
	if(l>=nx*ny) return 0;
	if((geti(l)==1)||(geti(l)==nx)||(getj(l)==1)||(getj(l)==ny)) return 0;
	if(l<0) return 0;
	return l;
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
	for(int l = 0; l< nx*ny; ++l) sumdiff+=pow(phi_0.at(l)-phi_1.at(l),2);
	return sqrt(sumdiff);
}