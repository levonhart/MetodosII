#include <iostream>
#include "pvc.hpp"

#define A 3
#define B 7
#define C 7
#define D 7
#define E 3
#define F 1

#define P 1.0*(A + B + C)
#define T 1.0*(D + E + F)


#define EPSILON 0.1e-3

using namespace std;
using namespace Eigen;


double f(double x){
	return -(P)/(T);
}

double coefDx(double x){
	return 1.0/x;
}

double coefU(double x){
	return 0.0;
}

int main(int argc, char ** argv){
	double xi = 0.2, xf = 0.5,
			uContorno = 0.0;
	int nparticoes = 4 + (A + B + C + D + E + F)%4;
	VectorXd u(nparticoes-1);


	u = Contorno1D( xi, xf, uContorno, &coefDx, &coefU, &f, nparticoes);
	cout << u << endl;
	
	return 0;
}
