#include "pvi.hpp"
#include <iostream>
#include <cmath>

using namespace std;
using namespace Eigen;

MN_amb_pvi::MN_amb_pvi(int n, double t, Eigen::VectorXd S, VectorXd (*F)(double t, VectorXd S)){
	dim = n;
	t0 = t;
	S0 = S;
	f = F;
}


VectorXd MN_amb_pvi::euler(double t, int npassos){
	if (npassos < 1) return VectorXd::Constant(dim,NAN);
	double dt = (t - t0)/npassos;
	VectorXd Sn(dim);
	Sn = S0;

	for(int i=0; i<npassos; i++){
		Sn = Sn + dt*f(t0 + i*dt, Sn);
	}
	
	return Sn;
}

VectorXd MN_amb_pvi::eulerMod(double t, int npassos){
	double dt = (t - t0)/npassos;
	VectorXd Sn(dim);
	Sn = S0;

	for(int i=0; i<npassos; i++){
		Sn = Sn + dt*f(t0 + i*dt + 0.5*dt, Sn + 0.5*dt*f(t0+i*dt,Sn));
	}
	
	return Sn;
}


VectorXd MN_amb_pvi::eulerApr(double t, int npassos){
	double dt = (t - t0)/npassos;
	VectorXd Sn(dim);
	Sn = S0;

	for(int i=0; i<npassos; i++){
		Sn = Sn + 0.5*dt*(f(t0 + i*dt, Sn) + f(t0 + (i+1)*dt, Sn + dt*f(t0 + i*dt, Sn)));
	}
	
	return Sn;
}



VectorXd MN_amb_pvi::rungeKutta(double t, int npassos){
	double dt = (t-t0)/npassos;
	VectorXd Sn(dim), k0(dim), k1(dim), k2(dim), k3(dim);
	Sn = S0;

	for(int i=0; i<npassos; i++){
		k0 = dt*f(t0 + i*dt, Sn);
		k1 = dt*f(t0 + i*dt + 0.5*dt, Sn + 0.5*k0);
		k2 = dt*f(t0 + i*dt + 0.5*dt, Sn + 0.5*k1);
		k3 = dt*f(t0 + i*dt + dt, Sn + k2);

		Sn = Sn + (1.0/6)*(k0 + 2*k1 + 2*k2 + k3);
	}

	return Sn;
}


VectorXd MN_amb_pvi::eulerIpct(double t, int npassos){
	double dt = (t - t0)/npassos;
	int i,j;
	VectorXd Sn(dim), a(dim), b(dim), c(dim), h(dim),
					Fa(dim), Fb(dim), Fc(dim);
	MatrixXd J(dim,dim);
	Sn = S0;

	for(i=0; i<npassos; i++){
		b = Sn;
		h.setConstant(HUGE_VAL);
		while(h.norm() > 0.1*dt){  // metodo de newton para encontrar raizes de funcoes
			Fb = b - Sn - dt*f(t0 + (i+1)*dt, b);
			for (j = 0; j < dim; j++){
				c = a = b;
				a(j) -= 0.5*dt;
				Fa = a - Sn - dt*f(t0 + (i+1)*dt, a);

				c(j) += 0.5*dt;
				Fc = c - Sn - dt*f(t0 + (i+1)*dt, c);

				J.col(j) = (1/dt)*( Fc - Fa );
			}
			// cout << " Jacobiana\n" << J; 
			h = J.inverse()*Fb;
			b = b - h;
			// cout << "\nb " << b.transpose() << "\t\terro "<< Fb.norm() << "\n";
		}

		Sn = b;

	}
	
	return Sn;
}