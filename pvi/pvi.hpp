#ifndef MN_PVI
#define MN_PVI

#include <Eigen/Dense>

struct MN_amb_pvi{
  private:
	int dim;
	double t0;
	Eigen::VectorXd S0;
	Eigen::VectorXd (*f)(double t, Eigen::VectorXd S);


  public:
	MN_amb_pvi(int n, double t, Eigen::VectorXd S, Eigen::VectorXd (*f)(double t, Eigen::VectorXd S));
	Eigen::VectorXd euler(double t, int npassos);
	Eigen::VectorXd eulerMod(double t, int npassos); // euler modificado y1 = y0 + h f(   x + 0.5 h  ,   y0 + 0.5 h f(x0, y0)  )
	Eigen::VectorXd eulerApr(double t, int npassos); // euler aprimorado y1 = y0 + 0.5 h[ f(x , y0) + f(   x+h  ,   y0 + h f(x0, y0)  ) ]
	Eigen::VectorXd rungeKutta(double t, int npassos);

	Eigen::VectorXd eulerIpct(double t, int npassos); // backwards euler ou euler implicito

};



#endif
