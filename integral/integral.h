#ifndef MN_INTEGRAL
#define MN_INTEGRAL

#define MN_NWTC_ABERTO 0xA
#define MN_NWTC_FECHADO 0xF

struct MN_amb_integral{
  private:
  	double *param;
	double (*f)(double x, double* param);

	double mudVarLegendre(double x, double a, double b);

  public:
	MN_amb_integral(double (*f)(double x, double* param), double * param);
	double calc(double x);

	double newtonCotes(double a, double b, int grau, int nparticoes, char filosofia); // filosofia : MN_ABERTO ou MN_FECHADO
	double gaussLegendre(double a, double b, int npontos, int nparticoes);

	double gaussLaguerre(double a,int npontos); // funcoes da forma exp{-x}*f(x) de a até +inf
	double gaussHermite(int npontos); // funcoes da forma exp{-x²}*f(x) de -inf a +inf
	double gaussChebyshev(int npontos); // funcoes da forma f(x)/sqrt(1-x²)   de -1 a 1

	double exponencialSimples(int npontos, int nparticoes);
	double exponencialDupla(int npontos, int nparticoes);




};



#endif
