#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <iomanip>

#include "Plot.h"

using namespace std;

vector<double> getgrid(int n, double h, double x0);
vector<double> f(vector<double> x);
double rectangle(vector<double> f_x, double h);
double trapezi(vector<double> f_x, double h);
double simpson(vector<double> f_x, double h);

int main(){
	
	double x0=0;
	double h=0.1;
	double n[]={10, 100, 1000, 1000000};
	
	// open writing file
	ofstream file("simpson.txt");
	
	Plot pl;
	vector<double> grid, func;
	for(int i=0; i<4; ++i){
		grid = getgrid(n[i], h, x0);
		func = f(grid);
		pl.x.push_back(n[i]);
		pl.y.push_back(simpson(func,h));
		
		grid.clear();
		func.clear();
	}
	
	// writing data on file
	double diff;
	for(int i= 0; i < pl.x.size(); ++i){
		//diff=pl.y.at(i)-sqrt(M_PI)*0.5;
		file << pl.x.at(i) << '\t' << setprecision(20)<< pl.y.at(i) << endl;
	}
	return 0;
}