#ifndef MN_PVC
#define MN_PVC

#include <Eigen/Dense>

typedef double (*funcaoReal)(double x);
typedef double (*funcaoReal2D)(double x, double y);

Eigen::VectorXd Contorno1D(double xi, double xf, double uContorno, funcaoReal coefDx, 
					funcaoReal coefU, funcaoReal F, int nparticoes);

Eigen::VectorXd Contorno2D(double xi, double xf, double yi, double yf, double uContorno, double coefLapl,
					funcaoReal2D coefDxDy, funcaoReal2D coefDx, funcaoReal2D coefDy, 
					funcaoReal2D coefU, funcaoReal2D F, int nparticoesx, int nparticoesy);




#endif
