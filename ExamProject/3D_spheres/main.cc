#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>

using namespace std;

//POISSON EQUATION SOLVER spherical charges (3D)

double p0 = 10;             
int lc =   2;
int px =   5, py =   5, pz = 5;
int lx =  10, ly =  10, lz = 10;
int nx=  100, ny = 100, nz = 100; 
double hx;
double hy;
double hz;
bool debug = true;

double getMlm( int l, int m );
double getl( int i, int j, int k );
double getj( int l );
double geti( int l );
double getk( int l );
double computediff( vector<double> phi_0, vector<double> phi_1 );

int main(){
	
	hx = double(lx)/nx;
	hy = double(ly)/ny;
	hz = double(lz)/nz;
	
	/*ofstream file("p_r_test.txt");
	ofstream file1("phi_r_test.txt");
	
	
	
	// initializing charge distribution
	vector<double> p;
	p.reserve(nx*ny*nz);
	for(int l = 0; l< nx*ny*nz; ++l ) p.push_back(0);
	
	if( debug ) cout<<"charge distribution vector initialized"<<endl;
	
	// filling charge distribution
	int i,j,k;
	int l;
	for(int i=0; i<nx; ++i){
		for(int j=0; j<ny; ++j){
			for(int k=0; k<nz; ++k){
				l = getl(i,j,k);
				if((sqrt(pow(i-px*0.5/hx,2)+pow(j-py*0.5/hy,2)+pow(k-pz*0.5/hz,2))< (lc*0.5/hx+0.07))&&
				   (sqrt(pow(i-px*0.5/hx,2)+pow(j-py*0.5/hy,2)+pow(k-pz*0.5/hz,2))> (lc*0.5/hx-0.07))){
					p.at(l) = p0;
					file<<i<<' '<<j<<' '<<k<<endl;
				}else if((sqrt(pow(i-px*1.5/hx,2)+pow(j-py*1.5/hy,2)+pow(k-pz*1.5/hz,2))< (lc*0.5/hx+0.07))&&
				         (sqrt(pow(i-px*1.5/hx,2)+pow(j-py*1.5/hy,2)+pow(k-pz*1.5/hz,2))> (lc*0.5/hx-0.07)))  {
					p.at(l) = p0;
					file<<i<<' '<<j<<' '<<k<<endl;
				}
			}
		}
		file<<endl;
	}
	
	if( debug ) cout<<"charge distribution vector filled"<<endl;
	
	//double tolerance = 0.1;
	// first iteration
	vector<double> phi_0, phi_1;
	phi_0.reserve(nx*ny*nz);
	phi_1.reserve(nx*ny*nz);
	for(int l = 0; l< nx*ny*nz; ++l ){
		phi_0.push_back(-p.at(l)/getMlm(l,l));
		phi_1.push_back(-p.at(l)/getMlm(l,l));
	} 
	if( debug ) cout<<"potential distribution vector initialized correctly"<<endl;
	
	// iterate
	int n= 0;
	double f1,fn,diff;
	f1= pow(hx,-2);
	fn= -2.*(pow(hx,-2)+pow(hy,-2)+pow(hz,-2));
	do{
		cout<<n<<" ";
		for(int l = 0; l< nx*ny*nz; ++l ){
			phi_0.at(l)=phi_1.at(l);
			phi_1.at(l)=0;
		} 
		for(int m = 0; m< nx*ny*nz; ++m ){
			i=geti(m);
			j=getj(m);
			k=getk(m);
			phi_1.at(m) = -f1*phi_0.at(getl(i+1,j,k))/fn;
			phi_1.at(m) += -f1*phi_0.at(getl(i-1,j,k))/fn;
			phi_1.at(m) += -f1*phi_0.at(getl(i,j+1,k))/fn;
			phi_1.at(m) += -f1*phi_0.at(getl(i,j-1,k))/fn;
			phi_1.at(m) += -f1*phi_0.at(getl(i,j,k+1))/fn;
			phi_1.at(m) += -f1*phi_0.at(getl(i,j,k-1))/fn;
			phi_1.at(m) += -p.at(m)/fn;
	    }
		diff=computediff(phi_0, phi_1);
		cout<<diff<<endl;
		n++;
	}while(n<100);//diff>0.1);

	cout<<diff<<" "<<endl;
	for(int l = 0; l< nx*ny*nz; ++l ){
		file1<<geti(l)<<' '<<getj(l)<<' '<<getk(l)<<' '<<phi_1.at(l)<<endl;
	}*/
	
	ifstream filein("phi_r.txt");
	vector<double>  phi_1;
	double r;
	int i,j,k;
	phi_1.reserve(nx*ny*nz);
	for(int l = 0; l< nx*ny*nz; ++l ){
		filein>>r;
		filein>>r;
		filein>>r;
		filein>>r;
		phi_1.push_back(r);
	}
	
	//electric field
	ofstream file2("Efield_test.txt");
	double ex,ey,ez,mod;
	for(int l = 0; l< nx*ny*nz; ++l ){
		i=geti(l);
		j=getj(l);
		k=getk(l);
		if((geti(l)==1)||(geti(l)==nx)||(getj(l)==1)||(getj(l)==ny)||(getk(l)==1)||(getk(l)==nz)) continue;
		
		ex = -(phi_1.at(getl(i+1,j,k))-phi_1.at(getl(i-1,j,k)))/hx;
		ey = -(phi_1.at(getl(i,j+1,k))-phi_1.at(getl(i,j-1,k)))/hy;
		ez = -(phi_1.at(getl(i,j,k+1))-phi_1.at(getl(i,j,k-1)))/hz;
		mod = sqrt(pow(ex,2)+pow(ey,2)+pow(ez,2))*0.2;
		ex /= 1.*mod;
		ey /= 1.*mod;
		ez /= 1.*mod;
	
		file2<<i<<' '<<j<<' '<<k<<' '<<ex<<' '<<ey<<' '<<ez<<endl;
		
		bool flag = false;
		for(int i =0; i<14; ++i){
			if(l%(nx*ny-1)==0){
				l+= 14*nx*ny;
				flag = true;
				break;
			}
			if(l%(nx-1) == 0){
				l+= 14*nx;
				flag = true;
				break;
			} 
			l++;
		}
		if(flag) continue;
		if(l%(nx*ny-1)==0) l+= 14*nx*ny;
		if(l%(nx-1) == 0) l+= 14*nx;
	}
	return 0;
}