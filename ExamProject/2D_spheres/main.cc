#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>

using namespace std;

//POISSON EQUATION SOLVER FOR spherical charges distribution (2D)

double p0 = 10;             // charge value
int lc =   2;				// diameter of spherical distribution
int px =   5, py =   5;     // middle point between the spheres
int lx =  10, ly =  10;     // dimensions of the system
int nx=  100, ny = 100;     // number of point for each dimension
double hx;
double hy;
bool debug = true;

double getMlm( int l, int m );                                      // get the element (l,m) of the matrix used for Jacobi method
double getl( int i, int j );                                        // get l index from i,j
double getj( int l );												// get j index from l
double geti( int l );												// get i index from l
double computediff( vector<double> phi_0, vector<double> phi_1 );   // compute difference between iterations

int main(){
	
	hx = double(lx)/nx;   // spacing of x grid
	hy = double(ly)/ny;   // spacing of y grid
	
	ofstream file("p_r_test.txt");    // charge distribution output file
	ofstream file1("phi_r_test.txt"); // potential output file
	
	
	
	// initializing charge distribution
	vector<double> p;
	p.reserve(nx*ny);
	for(int l = 0; l< nx*ny; ++l ) p.push_back(0);
	
	if( debug ) cout<<"charge distribution vector initialized"<<endl;
	
	// filling charge distribution
	int i,j;
	for(int l = 0; l< nx*ny; ++l  ){
		i = geti(l);
		j = getj(l);
		
		if(sqrt(pow(i-px*0.5/hx,2)+(pow(j-py*0.5/hy,2)))<lc*0.5/hx) p.at(l) = p0;
		else if(sqrt(pow(i-(px*1.5)/hx,2)+(pow(j-py*1.5/hy,2)))<lc*0.5/hx) p.at(l) = -p0;
		file<<i<<' '<<j<<' '<<p.at(l)<<endl;
	}	
	if( debug ) cout<<"charge distribution vector filled"<<endl;

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
	f1= pow(hx,-2);
	fn= -2.*(pow(hx,-2)+pow(hy,-2)); // note that fn == M_ll
	do{
		cout<<n<<" ";
		for(int l = 0; l< nx*ny; ++l ){
			phi_0.at(l)=phi_1.at(l);
			phi_1.at(l)=0;
		} 
		/*for(int m = 0; m< nx*ny; ++m ){
			i=geti(m);
			j=getj(m);
			phi_1.at(m) = -f1*phi_0.at(getl(i+1,j))/fn;
			phi_1.at(m) += -f1*phi_0.at(getl(i-1,j))/fn;
			phi_1.at(m) += -f1*phi_0.at(getl(i,j+1))/fn;
			phi_1.at(m) += -f1*phi_0.at(getl(i,j-1))/fn;
			phi_1.at(m) += -p.at(m)/fn;
	    }*/
		for(int m = 0; m< nx*ny; ++m ){
			i=geti(m);
			j=getj(m);
			phi_1.at(m) = -f1*phi_0.at(getl(i+1,j))/getMlm(m,m);
			phi_1.at(m) += -f1*phi_0.at(getl(i-1,j))/getMlm(m,m);
			phi_1.at(m) += -f1*phi_0.at(getl(i,j+1))/getMlm(m,m);
			phi_1.at(m) += -f1*phi_0.at(getl(i,j-1))/getMlm(m,m);
			phi_1.at(m) += -p.at(m)/getMlm(m,m);
	    }
		diff=computediff(phi_0, phi_1);
		cout<<diff<<endl;
		n++;
	}while(diff>0.0001);

	cout<<diff<<" "<<endl;;
	for(int l = 0; l< nx*ny; ++l ){
		file1<<geti(l)<<' '<<getj(l)<<' '<<phi_1.at(l)<<endl;
	}
	/*ifstream filein("phi_r.txt");
	vector<double>  phi_1;
	double r;
	int i,j;
	phi_1.reserve(nx*ny);
	for(int l = 0; l< nx*ny; ++l ){
		filein>>r;
		filein>>r;
		filein>>r;
		phi_1.push_back(r);
	}*/
	
	// compute electric field
	ofstream file2("Efield_test.txt");
	double ex,ey,mod;
	
	for(int i=0; i<nx; i++){
		for(int j=0; j<ny; j++){
			ex = -(phi_1.at(getl(i+1,j))-phi_1.at(getl(i-1,j)))/hx;
			ey = -(phi_1.at(getl(i,j+1))-phi_1.at(getl(i,j-1)))/hy;
			file2<<i<<' '<<j<<' '<<phi_1[getl(i,j)]<<' '<<ex<<' '<<ey<<endl;
		}
		file2<<endl;
	}
	return 0;
}