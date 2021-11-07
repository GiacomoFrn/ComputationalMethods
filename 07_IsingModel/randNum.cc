#include <iostream>
#include <cmath>

double randNum(){
	
    const long int   a = 16807;
    const long int   c = 0;
    const long       m = 2147483647;
    static long int x0 = 1;
    long int x1;
    double    r;
    
    x1=(a*x0+c)%m;
    r=((double) x1)/((double) m);
    x0=x1;
    
    return r;
    
}