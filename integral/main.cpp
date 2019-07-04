#include <iostream>
#include <iomanip>
#include "integral.h"
#include <cmath>

double fun(double x , double *param){
	return 4*x*x;
}

int main(){
	MN_amb_integral amb(fun,NULL);
	
	int nparticoes = 5;
	double intervalo_inferior_de_integracao = 0.0,
		intervalo_superior_de_integracao = 4.0;

	// Newton-cotes
	int grau_newtonCotes = 6;
	char tipoNewtonCotes = MN_NWTC_FECHADO; // MN_NWTC_FECHADO ou MN_NWTC_ABERTO
	
	// Gauss-Legendre
	int npontosGaussLegendre = 4;

	// Gauss-Laguerre / Gauss-Hermite / Gauss-Chebyshev
	int npontos_especiais = 4;

	std::cout << std::setprecision(30);
	std::cout << "neton-cotes:\t" << amb.newtonCotes(intervalo_inferior_de_integracao,
								intervalo_superior_de_integracao,
								grau_newtonCotes,
								nparticoes,
								tipoNewtonCotes) << std::endl;
	std::cout <<  "gauss-legendre:\t" << amb.gaussLegendre(intervalo_inferior_de_integracao,
								intervalo_superior_de_integracao,
								npontosGaussLegendre,
								nparticoes) << "\n\n\n";


	std::cout <<  "gauss-laguerre:\t" << amb.gaussLaguerre(0, npontos_especiais) << std::endl;
	std::cout <<  "gauss-hermite:\t" << amb.gaussHermite(npontos_especiais) << std::endl;
	std::cout <<  "gauss-chebyshev:\t" << amb.gaussChebyshev(npontos_especiais) << std::endl;

	// intervalo de integracao (a,b)
	// newtonCotes( a,  b, grau, nº de particoes, filosofia); // filosofia : MN_NWTC_ABERTO ( 0xA )  ou MN_NWTC_FECHADO  ( 0xF )
	// gaussLegendre( a, b, nº de pontos, nº de particoes);

	// gaussLaguerre( a, nº de pontos); 						// funcoes da forma exp{-x}*f(x) de a até +inf
	// gaussHermite( nº de pontos ); 						// funcoes da forma exp{-x²}*f(x) de -inf a +inf
	// gaussChebyshev( nº de pontos ); 						// funcoes da forma f(x)/sqrt(1-x²)   de -1 a 1

	// exponencialSimples( nº de pontos, nº de particoes);
	// exponencialDupla( nº de pontos, nº de particoes);
}
