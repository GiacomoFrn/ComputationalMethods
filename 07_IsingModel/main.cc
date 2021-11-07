#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>

#include "IsingLattice.h"

using namespace std;

double randNum();

int main(){
	
	int     Nx = 50;
	int     Ny = 50;
	double   J = 1;
	int Nsteps = 1000000;
	int kbT    = 1;
	
	vector<double> en, magn;
	magn.reserve(Nsteps);
	en.reserve(Nsteps);
	
	double ene;
	
	ofstream file("spinMat1.txt");
	
	IsingLattice* lattice2D = new IsingLattice( Nx, Ny, J );
	lattice2D->totalShuffle();
	
	ene = lattice2D->getEn();
	en.push_back(ene);
	magn.push_back(lattice2D->getMagn());
	
	int i = 0;
	double diffen;
	double probRatio;
	double r,r1;
	
	while( i<Nsteps ){
		
		if(i%1000==0) cout<<"passo: "<<i<<endl;
		
		r = randNum();
		lattice2D->singleShuffle( r );
		diffen = lattice2D->getDiffEn();
		probRatio = exp( -diffen/kbT );
		// Metropolis algorithm
		if( probRatio>1. ){
			ene += diffen;
			en.push_back(ene);
			magn.push_back(lattice2D->getMagn());
		}else{
			r1 = randNum();
			if( r1 < probRatio ){
				ene += diffen;
				en.push_back(ene);
				magn.push_back(lattice2D->getMagn());
			}else{
				lattice2D->singleShuffle( r );
				en.push_back(ene);
				magn.push_back(lattice2D->getMagn());
			}
		}
		
		if(i%100==0) lattice2D->writeSpinMat(file);
		if(i%100==0) file<<endl<<endl;
		i++;
	}
	
	ofstream out( "MC1.txt");
	out.precision(10);
	
	double avEn     = 0;
	double avMagn   = 0;
	double avSqEn   = 0;
	double avSqMagn = 0;
	
	for ( i=0; i<Nsteps; ++i){
		avEn     += en[i]/(Nx*Ny);
		avMagn   += magn[i]/(Nx*Ny);
        avSqEn   += pow(en[i]/(Nx*Ny),2);
        avSqMagn += pow(magn[i]/(Nx*Ny),2);
		if(i%100==0) out <<i<<" "<<en[i]/(Nx*Ny)<<" "<<magn[i]/(Nx*Ny)<<endl;
	}
	
	avEn /= Nsteps;
	avMagn /= Nsteps;
	avSqEn /= Nsteps;
	avSqMagn /= Nsteps;
	
	cout<< "Energia per sito: "<<avEn<<" +- "
	    <<sqrt((avSqEn-pow(avEn,2))/Nsteps)<<endl;
	cout<< "Magnetizzazione per sito: "<<avMagn<<" +- "
	    <<sqrt((avSqMagn-pow(avMagn,2))/Nsteps)<<endl;
	
	return 0;
}