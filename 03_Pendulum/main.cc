#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
using namespace std;
const double l= 1000;
const double g = 9.806/l;
struct data{
	
	double theta;
	double omega;
	
};
data func( data y){
	data tmp;
	tmp.theta = y.omega;
	tmp.omega = -g*sin(y.theta);
	return tmp;
}
data sc_prod (data y, double a){
	data tmp;
	tmp.omega = y.omega*a;
	tmp.theta = y.theta*a;
	return tmp;
}
data sum (data x, data y){
	data tmp;
	tmp.omega = y.omega+x.omega;
	tmp.theta = y.theta+x.theta;
	return tmp;
}
data RK( data y_i, double deltaT ){
	data Y_1 = y_i;
	data Y_2 = sum(y_i, sc_prod(func(Y_1), deltaT*0.5));
	data Y_3 = sum(y_i, sc_prod(func(Y_2), deltaT*0.5));
	data Y_4 = sum(y_i, sc_prod(func(Y_3), deltaT));
	return sum(y_i, sc_prod(sum(func(Y_1), sum(sc_prod(func(Y_2),2), sum(sc_prod(func(Y_3),2),func(Y_4)))), deltaT/6));
}
// positivo = true, negativo false
bool sign(double x){
	return (x>0);
}
int main(){
	double th_tau=2*M_PI/sqrt(g);
	double deltaT = 0.001;
	double Tmax = th_tau/2;
	int dim = Tmax/deltaT +1;
	//double t=0;
	//ofstream file1("omega.txt");
	//ofstream file1("theta.txt");
	ofstream file2("periodi_1000.txt");
	
	double m, c, tau;
	double theta_0=0.01;
	while(theta_0<(M_PI*0.5)){
		vector<data> Y(dim);
		data tmp,tmp1;
		tmp.omega = 0;
		tmp.theta = theta_0;
		Y.at(0)= tmp;
		//file1<< "0" <<'\t'<< Y.at(0).theta << endl;
		for(unsigned int i = 1; i<Y.size(); i++){
			tmp = RK( Y.at(i-1), deltaT);
			//file1<< deltaT*i <<'\t'<< tmp.theta << endl;
			//file1<< t+deltaT*i <<'\t'<< tmp.omega << endl;
			Y.at(i) = tmp;
		}
		for(unsigned int i = 1; i<Y.size(); i++){
			tmp1=Y.at(i);
			tmp= Y.at(i-1);
			if(sign(tmp.theta)!= sign(tmp1.theta)){
				m= (tmp1.theta-tmp.theta)/deltaT;
				c= tmp1.theta - m*(deltaT*i);
				break;
			}
		}
		tau=-c*4/m;
		file2<<theta_0<<'\t'<<tau/th_tau<<endl;
		Y.clear();
		theta_0+=0.01;
	}
	return 0;
	
	
}