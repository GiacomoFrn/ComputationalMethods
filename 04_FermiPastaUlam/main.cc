#include <iostream>
#include <fstream>
#include <vector>

#include "System.h"

using namespace std;

int main(){
	
	int         N = 31; // number of oscillators without the two fixed extremals
	double      k = 1; // hooke constant
	double      m = 1; // mass of oscillators
    double deltaT = 0.8; // delta time
	double      A = 0.1; // amplitude
	double      a = 1; // quadratic coefficient
	double 		b = 0; // 3^ coefficient
	
	System s(N, k, m, deltaT, A, a, b);
	ofstream file1("En_1.txt");
	ofstream file2("En_2.txt");
	ofstream file3("En_3.txt");
	ofstream file4("En_4.txt");
	ofstream file5("En_5.txt");
	
	//ofstream pos("Pos_gif.txt");

	
	double energy;
	double t=0;
	double Ncycles = 0;
	
	while( t < 1000000000000000 ){
		/*if(10*Ncycles<=t){
			for(int i=0; i<N+2;++i){
				pos<<i<<" "<<s.getPrPos(i)<<endl;
			}
			pos<<endl<<endl;
			Ncycles++;
		}*/
	
		if(100000*Ncycles<=t){
			
		
		s.computeEnergy(1);                     // calcolo energia modo n
		energy = s.getEnergy(1);    		// salvo energia modo n
		file1 << t << ' ' << energy << endl;
		s.computeEnergy(2);
		energy = s.getEnergy(2);
		file2 << t << ' ' << energy << endl;
		s.computeEnergy(3);
		energy = s.getEnergy(3);
		file3 << t << ' ' << energy << endl;
		s.computeEnergy(4);
		energy = s.getEnergy(4);
		file4 << t << ' ' << energy << endl;
		s.computeEnergy(5);
		energy = s.getEnergy(5);
		file5 << t << ' ' << energy << endl;
		
		Ncycles++;
		
		}
		
		//cout<< s.getEnergy(1) << ' ' << s.getEnergy(2)<< endl;    // stampo energia modo n
		//cout<<t<<" "<<s.getPrPos(1)<<endl;
		t+=deltaT;                            // incremento tempo
		
		s.Compute();      // calcolo pos->ne e vel->ne
		s.Update();       // salvo pos/vel -> ne in pr
	}
	
	
	
	return 0;
	
}