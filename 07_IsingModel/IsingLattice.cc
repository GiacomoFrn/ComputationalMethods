#include "IsingLattice.h"

#include <cmath>
#include <iostream>

using namespace std;

double randNum();

IsingLattice::IsingLattice(int nx, int ny, double j): 
	// inizializations
Nx(nx), Ny(ny), J(j){
	
	// inizial configuration of spins
	spinMat.reserve(Nx*Ny);
	for( int i=0; i<Nx*Ny; ++i){
		spinMat.push_back(1);
	}
}


IsingLattice::~IsingLattice() {
}

// shuffle a single spin in lattice and save diffEn
// respect to previous configuration
void IsingLattice::singleShuffle( double r ){
	
	int    l = floor(r*Nx*Ny);  // l position spin to be inverted
	int    j = floor(l*1./Nx);
	int    i = l-j*Nx;
	//cout<<i<<" "<<j<<endl;
	int nearX, nearY, nearL;
	double en = 0;
	
	// right X
	nearX = ( i<Nx-1 )? i+1 : 0; // periodic boundary conditions
	nearL = j*Nx + nearX;
	en += -J*(-1)*spinMat[l]*spinMat[nearL]; // H trial
	en -= -J*spinMat[l]*spinMat[nearL];      // H0
	
	// left X
	nearX = ( i>0 )? i-1 : Nx-1; // periodic boundary conditions
	nearL = j*Nx + nearX;
	en += -J*(-1)*spinMat[l]*spinMat[nearL]; // H trial
	en -= -J*spinMat[l]*spinMat[nearL];      // H0
	
	// up Y
	nearY = ( j<Ny-1 )? j+1 : 0; // periodic boundary conditions
	nearL = nearY*Nx + i;
	en += -J*(-1)*spinMat[l]*spinMat[nearL]; // H trial
	en -= -J*spinMat[l]*spinMat[nearL];      // H0
	
	// down Y
	nearY = ( j>0 )? j-1 : Ny-1; // periodic boundary conditions
	nearL = nearY*Nx + i;
	en += -J*(-1)*spinMat[l]*spinMat[nearL]; // H trial
	en -= -J*spinMat[l]*spinMat[nearL];      // H0
	
	spinMat[l] *= -1;
	
	diffEn = en;

	return;
	
}

void IsingLattice::totalShuffle(){
	
	double r;
	int    l = 0;
	
	for( int j=0; j<Ny; ++j){
		for( int i=0; i<Nx; ++i){
			r = randNum();
			spinMat[l] = ( r<0.5 )? -1 : +1;
			l++;
		}
	}
	
	return;
}

void IsingLattice::printSpinMat(){
	
	for( int j=0; j<Ny; ++j){
		for( int i=0; i<Nx; ++i){
			if ( i!=0) cout<<' ';
			cout<<spinMat[j*Nx+i];
		}
		cout<<endl;
	}
	
	return;
}

void IsingLattice::writeSpinMat(std::ofstream &file){
	
	for( int j=0; j<Ny; ++j){
		for( int i=0; i<Nx; ++i){
			if ( i==0) file<<spinMat[j*Nx+i];
			file<<' '<<spinMat[j*Nx+i];
		}
		file<<endl;
	}
}

double IsingLattice::getEn() const{
	
	double en = 0;
	int nearX, nearY, nearL;
	int l = 0;
	
	for( int j=0; j<Ny; ++j){
		for( int i=0; i<Nx; ++i){
	
			// right X
			nearX = ( i<Nx-1 )? i+1 : 0; // periodic boundary conditions
			nearL = j*Nx + nearX;
			en += -0.5*J*spinMat[l]*spinMat[nearL]; 
	
			// left X
			nearX = ( i>0 )? i-1 : Nx-1; // periodic boundary conditions
			nearL = j*Nx + nearX;
			en += -0.5*J*spinMat[l]*spinMat[nearL]; 
	
			// up Y
			nearY = ( j<Ny-1 )? j+1 : 0; // periodic boundary conditions
			nearL = nearY*Nx + i;
			en += -0.5*J*spinMat[l]*spinMat[nearL]; 
	
			// down Y
			nearY = ( j>0 )? j-1 : Ny-1; // periodic boundary conditions
			nearL = nearY*Nx + i;
			en += -0.5*J*spinMat[l]*spinMat[nearL];
			
			l++;
		}
	}
	
	return en;
}

double IsingLattice::getDiffEn() const{
	
	return diffEn;
}

double IsingLattice::getMagn() const{
	
	double mSum = 0;
	for( int m : spinMat ){
		mSum += m;
	}
	
	return mSum;
}
