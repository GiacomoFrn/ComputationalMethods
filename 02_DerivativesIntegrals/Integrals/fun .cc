#include <iostream>
#include <vector>
#include <cmath>
using namespace std;
// create equispaced grid given x0, h, n
vector<double> getgrid(int n, double h, double x0){
	vector<double> vec;
	for(int i=0;i<n;i++){
		vec.push_back(x0+i*h);
	}
	return vec;
}

// calculate function
vector<double> f(vector<double> x){
	vector<double> vec;
	for(int i=0;i<x.size();i++){
		vec.push_back(exp(-pow(x.at(i),2)));
	}
	return vec;	
}

// triangoli naif
double rectangle(vector<double> f_x, double h){
	double result=0;
	for(int i=0;i<f_x.size();i++){
		result += f_x.at(i)*h;
	}
	return result;
}

//trapezi
double trapezi(vector<double> f_x, double h){
	double result=0;
	for(int i=0;i<f_x.size()-1;i++){
		result += (f_x.at(i)+f_x.at(i+1))*h*0.5;
	}
	return result;
}

//simpson
double simpson(vector<double> f_x, double h){
	double result=0;
	for(int i=1;i<f_x.size()-1;i=i+2){
		result += ((1/3)*f_x.at(i-1)+(4/3)*f_x.at(i)+(1/3)*f_x.at(i+1))*2*h;
	}
	return result;
}