#include "integral.h"
#include <cmath>
#include <iostream>

MN_amb_integral::MN_amb_integral(double (*f)(double x, double *param), double * param){
	this->f = f;
	this->param = param;
}

double MN_amb_integral::mudVarLegendre(double s, double a, double b){
	return 0.5*( a+b + s*(b-a) );
}

double MN_amb_integral::newtonCotes(double a, double b, int grau, int nparticoes, char filosofia){
	double 	I = 0,
			dx = (b-a)/nparticoes,
			h, frac, x, *w;
	int i,j;
	if (grau == 0){
		double pesos[] = {1.0};
		frac = 2.0;
		w = pesos;
	} else if (grau == 1){
		double pesos[] = {1.0, 1.0};
		frac = 0.5;
		w = pesos;
	} else if (grau == 2 && filosofia == MN_NWTC_FECHADO){
		double pesos[] = {1.0, 4.0, 1.0};
		frac = 1.0/3.0;
		w = pesos;
	} else if (grau == 2 && filosofia == MN_NWTC_ABERTO){
		double pesos[] = {2.0, -1.0, 2.0};
		frac = 4.0/3.0;
		w = pesos;
	} else if (grau == 3 && filosofia == MN_NWTC_FECHADO){
		double pesos[] = {1.0, 3.0, 3.0, 1.0};
		frac = 3.0/8.0;
		w = pesos;
	} else if (grau == 3 && filosofia == MN_NWTC_ABERTO){
		double pesos[] = {11.0, 1.0, 1.0, 11.0};
		frac = 5.0/24.0;
		w = pesos;
	} else if (grau == 4 && filosofia == MN_NWTC_FECHADO){
		double pesos[] = {7.0, 32.0, 12.0, 32.0, 7.0};
		frac = 2.0/45.0;
		w = pesos;
	} else if (grau == 4 && filosofia == MN_NWTC_ABERTO){
		double pesos[] = {11.0, -14.0, 26.0, -14.0, 11.0};
		frac = 3.0/10.0;
		w = pesos;
	}

	 else if (grau == 5 && filosofia == MN_NWTC_FECHADO){
		double pesos[] = {19.0, 75.0, 50.0, 50.0, 75.0, 19.0};
		frac = 5.0/288.0;
		w = pesos;
	} else if (grau == 5 && filosofia == MN_NWTC_ABERTO){
		double pesos[] = {611.0, -453.0, 562.0, 562.0, -453.0, 611.0};
		frac = 7.0/1440.0;
		w = pesos;
	} else if (grau == 6 && filosofia == MN_NWTC_FECHADO){
		double pesos[] = {41.0, 216.0, 27.0, 272.0, 27.0, 216.0, 41.0};
		frac = 1.0/140.0;
		w = pesos;
	} else if (grau == 6 && filosofia == MN_NWTC_ABERTO){
		double pesos[] = {460.0, -954.0, 2196.0, -2459.0, 2196.0, -954.0, 460.0};
		frac = 8.0/945.0;
		w = pesos;
	} else if (grau == 7 && filosofia == MN_NWTC_FECHADO){
		double pesos[] = {751.0, 3577.0, 1323.0, 2989.0, 1323.0, 3577.0, 751.0};
		frac = 7.0/17280.0;
		w = pesos;
	// } else if (grau == 7 && filosofia == MN_NWTC_ABERTO){
	// 	double w[] = {11.0, -14.0, 26.0, -14.0, 11.0};
	// 	frac = 3.0/10.0;
		// w = pesos;
	} else if (grau == 8 && filosofia == MN_NWTC_FECHADO){
		double pesos[] = {989.0, 5888.0, -928.0, 10496.0, -4540.0, 10496.0, -928.0, 5888.0, 989.0};
		frac = 4.0/14175.0;
		w = pesos;
	// } else if (grau == 8 && filosofia == MN_NWTC_ABERTO){
	// 	double w[] = {11.0, -14.0, 26.0, -14.0, 11.0};
	// 	frac = 3.0/10.0;
		// w = pesos;
	} else {return NAN;};
	

	if(filosofia == MN_NWTC_FECHADO){
		if(grau < 1) return NAN;
		h = dx/(grau);
		x = a;
	} else if(filosofia == MN_NWTC_ABERTO){
		h = dx/(grau+2);
		x = a+h;
	} else{
		return NAN;
	}

	for (i = 0; i < nparticoes; i++){
		for (j = 0; j < grau+1; ++j){
			I += w[j]*f(x+j*h , param);
		}
		x += dx;
	}
	I *= h*frac;

	return I;	

}

double fun(double x , double *param){
	return sin(x);
}

int main(){
	MN_amb_integral amb(fun,NULL);
	printf("%50.40lf\n", amb.newtonCotes(0.0,1.0,6 ,5,MN_NWTC_ABERTO) );
}