#include <iostream>
#include <vector>
#include <cmath>
#include <complex.h>

using namespace std;

extern int        L;
extern double    x0;
extern double     q;
extern double sigma;
extern double     a;
extern double     b;
extern double    v0;
extern int       nx;
extern double    dt;
extern double    nt;
extern double     h;
extern complex<double> imUn;

complex<double> f( int i ){
	
	double x = h*i;
	
	return exp(imUn*q*x)*exp(-0.5*pow(x-x0,2)/pow(sigma,2));
	
}

//trapezi
double trapezi( vector<complex<double>> f_x ){
	
	double result=0;
	vector<double> f;
	f.reserve(f_x.size());
	for(unsigned int i=0;i<f_x.size()-1;i++){
		f.push_back(pow(abs(f_x[i]),2));
	}
	for(unsigned int i=0;i<f.size()-1;i++){
		result += (f.at(i)+f.at(i+1))*h*0.5;
	}
	
	return result;
}

double Vi( int i){
	
	if( (i*h>b)||(i*h<a) ) return 0;
	return v0;
}

complex<double> Mii( int i ){
	
	return imUn*4.*h*h/dt-2.-2.*h*h*Vi(i);
	
}

complex<double> bi( int i, vector<complex<double>> phi0 ){
	
	return -phi0[i+1]+2.*phi0[i]-phi0[i-1]+imUn*4.*h*h*phi0[i]/dt
	       +2.*h*h*Vi(i)*phi0[i];
		   
}