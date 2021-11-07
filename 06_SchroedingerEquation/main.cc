#include <iostream>
#include <fstream>
#include <vector>
#include <complex.h>
#include <cmath>

using namespace std;

int        L = 500;
double    x0 = 200;
double     q = 2;
double sigma = 20;
double     a = 250;
double     b = 260;
double    v0 = 1.7;
int       nx = 50;
double    dt = 0.1;
double    nt = 2000;
double     h;
complex<double> imUn (0.0, 1.0);

complex<double> f( int i );
double trapezi( vector<complex<double>> f_x );
double Vi( int i );
complex<double> Mii( int i );
complex<double> bi( int i, vector<complex<double>> phi0 );

int main(){
	
	h = L*1./nx;
	
	vector<complex<double>> phi0, phi1;
	phi0.reserve(nx+1);
	phi1.reserve(nx+1);
	
	// inizializzo
	phi0.push_back(0);
	phi1.push_back(0);
	for( int i = 1; i < nx; ++i ){
		phi0.push_back(f(i));
		phi1.push_back(0);
	}
	phi0.push_back(0);
	phi1.push_back(0);
	
	// normalizzo e scrivo il primo passo
	ofstream file( "phi.txt" );
	double normFact = sqrt(trapezi( phi0 ));
	for( int i = 0; i < nx+1; ++i ){
		phi0[i] /= normFact;
		file<<i*h<<' '<<pow(abs(phi0[i]),2)<<endl;
	}
	file<<endl;
	
	vector<complex<double>> alpha, beta;
	alpha.reserve(nx+1);
	beta.reserve(nx+1);
	for( int i = 0; i<nx+1; ++i ){
		alpha.push_back(0);
		beta.push_back(0);
	}
	
	double t = 1;
	
	while( t<nt ){
		
		alpha[1] = -Mii(1);
		beta[1]  = bi(1, phi0);
	
		for( int i = 2; i<nx-1; ++i ){
			alpha[i] = -1./alpha[i-1] -Mii(i);
			beta[i] = bi(i, phi0)+beta[i-1]/alpha[i-1];
		}
		file<<nx*h<<' '<<pow(abs(phi1[nx]),2)<<endl;
		phi1[nx-1] = pow(1./alpha[nx-2] +Mii(nx-1), -1)*(bi(nx-1, phi0)+beta[nx-2]/alpha[nx-2]);
		file<<(nx-1)*h<<' '<<pow(abs(phi1[nx-1]),2)<<endl;
		for( int i = nx-2; i>0; --i){
			phi1[i] = phi1[i+1]/alpha[i] -beta[i]/alpha[i];
			file<<i*h<<' '<<pow(abs(phi1[i]),2)<<endl;
		}
		file<<0*h<<' '<<pow(abs(phi1[0]),2)<<endl;
		file<<endl<<endl;
		
		normFact = sqrt(trapezi( phi1 ));
		cout<<"passo"<<t<<" norma:"<<normFact<<endl;
		for( int i = 0; i<nx+1; ++i ){
		phi0[i]=phi1[i];
		}
		
		t++;
	}
	
	return 0;
	
}