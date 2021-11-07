#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>

using namespace std;

//POISSON EQUATION

double p0 = 80;             
int lc =   2;
int px =   4, py =   4;
int lx =  10, ly =  10;
int nx=  100, ny = 100; 
double hx;
double hy;
bool debug = true;

double getl( int i, int j );
bool inLcP( int i, int j);
bool inLcN( int i, int j);
double getj( int l );
double geti( int l );
double getMlm( int l, int m );
double computediff( vector<double> phi_0, vector<double> phi_1 );

int main(){
	ifstream filein("p_r2.txt");
	ifstream file1in("phi_2.txt");
	
	hx = double(lx)/nx;
	hy = double(ly)/ny;
	
	// initializing charge distribution
	vector<double> p;
	p.reserve(nx*ny);
	/*for(int l = 0; l< nx*ny; ++l ) p.push_back(0);
	
	if( debug ) cout<<"charge distribution vector initialized"<<endl;
	
	// filling charge distribution
	for( int i = 1; i <= nx; ++i ){
		for ( int j = 1; j <= ny; ++j ){
			if(inLcP(i,j)) p.at(getl(i,j)) = p0;
			if(inLcN(i,j)) p.at(getl(i,j)) = -p0;
			file<<i<<' '<<j<<' '<<p.at(getl(i,j))<<endl;
		}
	file<<endl;
	}
	if( debug ) cout<<"charge distribution vector filled"<<endl;*/
	
	double r;
	for(int l = 0; l< nx*ny; ++l ){
		filein>>r;
		filein>>r;
		filein>>r;
		p.push_back(r);
	}
	if( debug ) cout<<"charge distribution vector read correctly"<<endl;
	filein.close();
	
	// reading potential distribution vector
	vector<double> phi_0, phi_1;
	phi_0.reserve(nx*ny);
	phi_1.reserve(nx*ny);
	for(int l = 0; l< nx*ny; ++l ){
		file1in>>r;
		file1in>>r;
		file1in>>r;
		phi_0.push_back(r);
		phi_1.push_back(r);
	} 
	
	if( debug ) cout<<"potential distribution vector read correctly"<<endl;
	file1in.close();
	
	ofstream file1("phi_r2p.txt");
	
	/*//double tolerance = 0.1;
	// first iteration
	vector<double> phi_0, phi_1;
	phi_0.reserve(nx*ny);
	phi_1.reserve(nx*ny);
	for(int l = 0; l< nx*ny; ++l ){
		phi_0.push_back(-p.at(l)/getMlm(l,l));
		//file<<geti(l)<<' '<<getj(l)<<' '<<phi_0.at(l)<<endl;
		phi_1.push_back(-p.at(l)/getMlm(l,l));
	} */
	
	// iterate
	int n= 0;
	do{
		for(int l = 0; l< nx*ny; ++l ){
			phi_0.at(l)=phi_1.at(l);
			phi_1.at(l)=0;
		} 
		for(int l = 0; l< nx*ny; ++l ){
			for(int m = 0; m< nx*ny; ++m ){
				if( l == m ) continue;
				if(getMlm(l,m)==0) continue;
				if(phi_0.at(m)==0) continue;
				phi_1.at(l) += -getMlm(l,m)*phi_0.at(m)/getMlm(l,l);
				//if(getMlm(l,m)!=0) cout<<getMlm(l,m)<<endl;
			}
			phi_1.at(l) += -p.at(l)/getMlm(l,l);
	    }
		cout<<computediff(phi_0, phi_1)<<endl;
		n++;
	}while(n<10000);
	
	for(int l = 0; l< nx*ny; ++l ){
		file1<<geti(l)<<' '<<getj(l)<<' '<<phi_1.at(l)<<endl;
	} 
	
	return 0;
}