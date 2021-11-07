#include <iostream>
#include <cmath>
#include <vector>

using namespace std;

extern double p0;             
extern int lc;
extern int px;
extern int py;
extern int pz;
extern int lx;
extern int ly;
extern int lz;
extern int nx;
extern int ny;
extern int nz;
extern double hx;
extern double hy;
extern int hz;
extern bool debug;


// return k from l
double getk( int l ){
	return floor(l/(nx*ny))+1;
}

// return j from l
double getj( int l ){
	return floor((l-(getk(l)-1)*nx*ny)/nx)+1;
}

// return i from l
double geti( int l ){
	return l-(getk(l)-1)*nx*ny-(getj(l)-1)*nx+1;
}

// return l from i,j,k
double getl( int i, int j, int k ){
	//if( debug ) std::cout<<(j-1)*nx+i-1<<std::endl;
	double l=(k-1)*nx*ny+(j-1)*nx+i-1;
	if(l>=nx*ny*ny) return 0;
	if((geti(l)==1)||(geti(l)==nx)||(getj(l)==1)||(getj(l)==ny)||(getk(l)==1)||(getk(l)==nz)) return 0;
	if(l<0) return 0;
	return l;
}


// return matrix element (l,m)
double getMlm( int l, int m ){
	// check if is an edge point
	if((geti(l)==1)||(geti(l)==nx)||(getj(l)==1)||(getj(l)==ny)||(getk(l)==1)||(getk(l)==nz)){
		if(l==m) return 1;
		else return 0;
	}
	
	if(((geti(l)==geti(m)-1)&&(getj(l)==getj(m))&&(getk(l)==getk(m)))||((geti(l)==geti(m)+1)&&(getj(l)==getj(m))&&(getk(l)==getk(m)))) return pow(hx,-2);
	if(((getj(l)==getj(m)-1)&&(geti(l)==geti(m))&&(getk(l)==getk(m)))||((getj(l)==getj(m)+1)&&(geti(l)==geti(m))&&(getk(l)==getk(m)))) return pow(hy,-2);
	if(((getj(l)==getj(m))&&(geti(l)==geti(m))&&(getk(l)==getk(m)-1))||((getj(l)==getj(m))&&(geti(l)==geti(m))&&(getk(l)==getk(m)+1))) return pow(hz,-2);
	if(l==m) return -2.*(pow(hx,-2)+pow(hy,-2)+pow(hz,-2));
	return 0;
}

// compute difference between iterations
double computediff( vector<double> phi_0, vector<double> phi_1 ){
	double sumdiff = 0;
	for(int l = 0; l< nx*ny*nz; ++l) sumdiff+= abs(phi_0.at(l)-phi_1.at(l))/(nx*ny*nz);
	return sumdiff;
}