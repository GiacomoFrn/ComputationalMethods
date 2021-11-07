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

// calculate sin(vec)
vector<double> Sin(vector<double> x){
	vector<double> vec;
	for(int i=0;i<x.size();i++){
		vec.push_back(sin(x.at(i)));
	}
	return vec;	
}

// calculate cos(vec)
vector<double> Cos(vector<double> x){
	vector<double> vec;
	for(int i=0;i<x.size();i++){
		vec.push_back(cos(x.at(i)));
	}
	return vec;	
}

// d/dx O(h) dx e sx (N-1)
vector<double> fod_dx(vector<double> f_x, double h){
	vector<double> der;
	double der_value;
	for(int i=0;i<f_x.size()-1;i++){
		der_value= (f_x.at(i+1)-f_x.at(i))/h;
		der.push_back(der_value);
	}
	return der;	
}

vector<double> fod_sx(vector<double> f_x, double h){
	vector<double> der;
	double der_value;
	for(int i=1;i<f_x.size();i++){
		der_value= (f_x.at(i)-f_x.at(i-1))/h;
		der.push_back(der_value);
	}
	return der;	
}

// d/dx O(h^2) (N-2)
vector<double> fod_2(vector<double> f_x, double h){
	vector<double> der;
	double der_value;
	for(int i=1;i<f_x.size()-1;i++){
		der_value= (f_x.at(i+1)-f_x.at(i-1))/(2*h);
		der.push_back(der_value);
	}
	return der;	
}

// d/dx O(h^4) (N-4)
vector<double> fod_4(vector<double> f_x, double h){
	vector<double> der;
	double der_value;
	for(int i=2;i<f_x.size()-2;i++){
		der_value= (-f_x.at(i+2)+8*f_x.at(i+1)-8*f_x.at(i-1)+f_x.at(i-2))/(12*h);
		der.push_back(der_value);
	}
	return der;	
}
