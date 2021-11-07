#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
using namespace std;
vector<double> getgrid_a(int n, double h, double x0){
	vector<double> vec;
	for(int i=0;i<n;i++){
		vec.push_back(x0+i*h);
	}
	return vec;
}
vector<double> getgrid_b(double a, double b, int n){
	vector<double> vec;
	for(int i=0;i<=n;i++){
		vec.push_back(i*(b-a)/n+a);
	}
	return vec;
}
struct plot{
	vector<double> x;
	vector<double> y;
};
vector<double> Sin(vector<double> x){
	vector<double> vec;
	for(int i=0;i<x.size();i++){
		vec.push_back(sin(x.at(i)));
	}
	return vec;	
}
int main(){
	char c;
	int n;
	double h,x0,a,b;
	plot data;
	ofstream file("plot_data.txt");
	cout<<" a or b?"<<endl;
	cin>>c;
	if(c=='a'){
		cout<<"insert n, h, x0"<<endl;
		cin>>n>>h>>x0;
		data.x=getgrid_a(n,h,x0);
	}else{
		cout<<"insert a, b, n"<<endl;
		cin>>a>>b>>n;
		data.x=getgrid_b(a,b,n);
	}
	data.y=Sin(data.x);
	file<<"X"<<'\t'<<"Y"<<endl;
	for(int i=0;i<data.x.size();i++){
		file<<data.x.at(i)<<'\t'<<data.y.at(i)<<endl;
	}
}

