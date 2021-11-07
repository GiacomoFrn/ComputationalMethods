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
	
	
	hx = double(lx)/nx;
	hy = double(ly)/ny;
	/*
	ifstream filein("p_r1.txt");
	// initializing charge distribution
	vector<double> p;
	p.reserve(nx*ny);
	
	double r;
	for(int l = 0; l< nx*ny; ++l ){
		filein>>r;
		filein>>r;
		filein>>r;
		p.push_back(r);
	}
	if( debug ) cout<<"charge distribution vector read correctly"<<endl;
	filein.close();
	
	
	
	//double tolerance = 0.1;
	// first iteration
	vector<double> phi_0, phi_1;
	phi_0.reserve(nx*ny);
	phi_1.reserve(nx*ny);
	for(int l = 0; l< nx*ny; ++l ){
		phi_0.push_back(-p.at(l)/getMlm(l,l));
		phi_1.push_back(-p.at(l)/getMlm(l,l));
	} 
	if( debug ) cout<<"potential distribution vector initialized correctly"<<endl;
	// iterate
	int n= 0;
	double f1,fn,diff;
	int i,j;
	f1= pow(hx,-2);
	fn= -2.*(pow(hx,-2)+pow(hy,-2));
	do{
		//cout<<n<<" ";
		for(int l = 0; l< nx*ny; ++l ){
			phi_0.at(l)=phi_1.at(l);
			phi_1.at(l)=0;
		} 
		for(int m = 0; m< nx*ny; ++m ){
			i=geti(m);
			j=getj(m);
			phi_1.at(m) = -f1*phi_0.at(getl(i+1,j))/fn;
			phi_1.at(m) += -f1*phi_0.at(getl(i-1,j))/fn;
			phi_1.at(m) += -f1*phi_0.at(getl(i,j+1))/fn;
			phi_1.at(m) += -f1*phi_0.at(getl(i,j-1))/fn;
			phi_1.at(m) += -p.at(m)/fn;
	    }
		diff=computediff(phi_0, phi_1);
		//cout<<diff<<endl;
		n++;
	}while(diff>0.1);
	ofstream file1("phi_r12.txt");
	cout<<diff<<" "<<endl;;
	for(int l = 0; l< nx*ny; ++l ){
		file1<<geti(l)<<' '<<getj(l)<<' '<<phi_1.at(l)<<endl;
	} */
	ifstream filein("phi_r12.txt");
	vector<double>  phi_1;
	double r;
	int i,j;
	phi_1.reserve(nx*ny);
	for(int l = 0; l< nx*ny; ++l ){
		filein>>r;
		filein>>r;
		filein>>r;
		phi_1.push_back(r);
	}
	//electric field
	ofstream file2("Efield.txt");
	double ex,ey;
	for(int l = 0; l< nx*ny; ++l ){
		i = geti(l);
		j= getj(l);
		ex = (phi_1.at(getl(i+1,j))-phi_1.at(getl(i-1,j)))/hx;
		ey = (phi_1.at(getl(i,j+1))-phi_1.at(getl(i,j-1)))/hy;
		file2<<i<<' '<<j<<' '<<ex<<' '<<ey<<endl;
		if(l%nx == 0){
			l+= 2*nx;
			continue;
		} 
		l++;
		if(l%nx == 0){
			l+= 2*nx;
			continue;
		} 
		l++;
		if(l%nx == 0) l+= nx;
	}
	return 0;
}

	
	
	