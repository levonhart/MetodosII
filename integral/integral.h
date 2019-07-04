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
	double gaussLegendre();
	double gassHermite();
	double gaussLaguerre();
	double gaussChebyshev();
	double exponencialSimples();
	double exponencialDupla();




};



#endif