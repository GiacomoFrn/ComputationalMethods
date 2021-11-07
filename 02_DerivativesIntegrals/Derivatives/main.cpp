#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>

#include "Plot.h"

using namespace std;
vector<double> getgrid(int n, double h, double x0);
vector<double> Sin(vector<double> x);
vector<double> Cos(vector<double> x);
vector<double> fod_dx(vector<double> f_x, double h);
vector<double> fod_sx(vector<double> f_x, double h);
vector<double> fod_2(vector<double> f_x, double h);
vector<double> fod_4(vector<double> f_x, double h);

int main(){
	
	// open writing file
	ofstream file("plot_data_a.txt");
	
	// building grid
	int n=100;
	double x0=0,h=0.1;
	Plot* pl= new Plot;
	pl->x=getgrid(n, h, x0);
	
	// bulding sin and cos
	vector<double> sin= Sin(pl->x);
	vector<double> cos= Cos(pl->x);
	
	//setting offset grid variables
	int end=0;
	int beg=0;
	
	// choosing order of derivation
	char x;
	cout << " choose order of derivative:" << endl;
	cout << "a) O(h)_dx" << endl;
	cout << "b) O(h)_sx" << endl;
	cout << "c) O(h^2)" << endl;
	cout << "d) O(h^4)" << endl;
	cin>>x;
	switch(x){
		case 'a':
			pl->y=fod_dx(sin, h);
			end=1;
			break;
			
		case 'b':
			pl->y=fod_sx(sin, h);
			beg=1;
			break;
			
		case 'c':
			pl->y=fod_2(sin, h);
			beg=1;
			end=1;
			break;
			
		case 'd':
			pl->y=fod_4(sin, h);
			beg=2;
			end=2;
			break;
	}
	
	// writing data on file
	for(int i= beg; i < pl->x.size()-end; ++i){
		file << pl->x.at(i) << '\t' << pl->y.at(i-beg)-cos.at(i) << endl;
	}
	return 0;
}