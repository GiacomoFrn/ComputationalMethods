#include "System.h"
#include <iostream>
#include<cmath>

System::System(int Num, double kf, double mass, double dT, double Amp, double al, double be ):
// initializations
N(Num), k(kf), m(mass), deltaT(dT), A(Amp), a(al), b(be) {
	
	Npos.reserve(N+2);
	Nvel.reserve(N+2);
	nC.reserve(N+2);
	nDC.reserve(N+2);
	nE.reserve(N+2);
	
	w = sqrt((4*k)/m)*sin((M_PI)/(2*(N+1)));
	
	Initialize();
}


System::~System() {
}

// get n mode period
double System::getNperiod( int n ){
	return 2*M_PI/(sqrt(4*k/m)*sin(M_PI*n*0.5/(N+1)));
}
// initialize position and velocity of each oscillator at n=1 normal mode
void System::Initialize(){
	
	System::Position pPos;
	pPos.pr = 0;
	pPos.ne = 0;
	
	Npos.push_back(pPos);
	
	System::Velocity pVel;
	pVel.pr = 0;
	pVel.ne = 0;
	
	Nvel.push_back(pVel);
	
	nC.push_back(0);
	nDC.push_back(0);
	nE.push_back(0);

	int i;
	for( i=1; i <= N; ++i){
		
		pPos.pr = 0; // Asin(pi i / N + 1) sin ( - w *0) = 0 for all i-1
		pPos.ne = 0;
		Npos.push_back(pPos);
	
		pVel.pr = -w*A*sin((M_PI*i)/(N+1)); // * cos(0) = 1
		//cout<<-w*A*sin((M_PI*i)/(N+1))<<endl;
		pVel.ne = 0;
		Nvel.push_back(pVel);
		//cout<<Npos[i].pr<<" "<<Nvel[i].pr<<endl;
		//cout<<Npos[i]->pr<<endl;
		nC.push_back(0);
		nDC.push_back(0);
		nE.push_back(0);
	}
	pPos.pr = 0;
	pPos.ne = 0;
	Npos.push_back(pPos);
	
	pVel.pr = 0;
	pVel.ne = 0;
	Nvel.push_back(pVel);
	cout<<Nvel.size()<<endl;
	nC.push_back(0);
	nDC.push_back(0);
	nE.push_back(0);
	
	// print max and min periods 
	std::cout << getNperiod(1)  << std::endl;
	/*for (int i = 0; i<N+2; i++){
		cout<<Npos[i].pr<<" "<<Nvel[i].pr<<endl;
	}*/
	return;
	
}

// compute new positions and velocities with velocity verlet algorithm
void System::Compute() {
	//cout<<"computing new positions and velocities"<<endl;

	int i;
	for( i = 1; i <= N; ++i){
		Npos[i].ne = Npos[i].pr + Nvel[i].pr*deltaT + FunctionPr( i )*pow(deltaT,2)*0.5/m;
	}
	for( i = 1; i <= N; ++i){
		Nvel[i].ne = Nvel[i].pr + deltaT*0.5*( FunctionPr( i ) + FunctionNe( i ) )/m;
	//	cout<<deltaT*0.5*( FunctionPr( i ) + FunctionNe( i ) )/m<<endl;
	}
	//cout<<( FunctionPr( 1 ) + FunctionNe( 1 ) )/m<<endl;
	/*for( i = 1; i <= N; ++i){
		cout<<Npos[i].pr<<" "<<Nvel[i].pr<<endl;
	}*/
	return;
}

// update present data
void System::Update() {
	
	int i;
	
	for( i = 1; i <= N; ++i){
		Npos[i].pr = Npos[i].ne;
		Nvel[i].pr = Nvel[i].ne; 
	}

	return;
}

// compute F(t)
double System::FunctionPr( int i){
	return   k*(Npos[i+1].pr-2*Npos[i].pr+Npos[i-1].pr) 
		   + a*(pow(Npos[i+1].pr-Npos[i].pr,2)-pow(Npos[i].pr-Npos[i-1].pr,2))
		   + b*(pow(Npos[i+1].pr-Npos[i].pr,3)-pow(Npos[i].pr-Npos[i-1].pr,3));
}

// compute F(t+deltaT)
double System::FunctionNe( int i){
	return   k*(Npos[i+1].ne-2*Npos[i].ne+Npos[i-1].ne) 
		   + a*(pow(Npos[i+1].ne-Npos[i].ne,2)-pow(Npos[i].ne-Npos[i-1].ne,2))
		   + b*(pow(Npos[i+1].ne-Npos[i].ne,3)-pow(Npos[i].ne-Npos[i-1].ne,3));
}


// compute energiy
void System::computeEnergy( int n ){
	//cout<<"Computing energy"<<endl;
	int i; 
	nC[n]  = 0;
	nDC[n] = 0;
	nE[n]  = 0;
	for( i = 1; i <= N; ++i){
		nC[n] += (Npos[i].pr)*sin(M_PI*n*i/(N+1));
		//cout<<nDC[n]<<" "<<Nvel[i].pr<<" "<<sin(M_PI*n*i/(N+1))<<endl;
		nDC[n]+= (Nvel[i].pr)*sin(M_PI*n*i/(N+1));
	}
	nC[n] *= sqrt(2./(N+1));
	//cout<<nDC[n]<<" "<<sqrt(2./(N+1))<<endl;
	nDC[n]*= sqrt(2./(N+1));
	//cout<<nC[n]<<endl;
	//cout<<nDC[n]<<endl;
	nE[n] = 0.5*m*pow(nDC[n],2)+2*k*pow(nC[n]*sin(M_PI*n*0.5/(N+1)),2);
	//cout<<nE[n]<<endl;
	return;
}

// get energy
double System::getEnergy( int n ){
	return nE[n];
}

// get N, number of oscillators -2
double System::getN() const{
	return N;
}

double System::getPrVel( int i ) const{
	return Nvel[i].pr;
}


double System::getNeVel( int i ) const{
	return Nvel[i].ne;
}


double System::getPrPos( int i ) const{
	return Npos[i].pr;
}


double System::getNePos( int i ) const{
	return Nvel[i].ne;
}

