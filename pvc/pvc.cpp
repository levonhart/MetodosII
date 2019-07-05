#include "pvc.hpp"
#include <iostream>
#include <cmath>

using namespace std;
using namespace Eigen;

VectorXd Contorno1D(double xi, double xf, double uContorno, funcaoReal coefDx, 
					funcaoReal coefU, funcaoReal F, int nparticoes){
	int dim = nparticoes-1, i, j;
    double x, dx = (xf - xi)/nparticoes;
    VectorXd u(dim),
             b(dim);
    MatrixXd A(dim, dim);

    for(i=0; i<dim; i++){
		x = xi+dx + i*dx;
		b(i) = F(x);

		A(i,i) = -(2.0/(dx*dx));
        if(i-1 >= 0) A(i,i-1) = 1.0/(dx*dx) - coefDx(x)/(2*dx);
        if(i+1 < dim) A(i,i+1) = 1.0/(dx*dx) + coefDx(x)/(2*dx);
	}
	b(0) -= (1.0/(dx*dx) - coefDx(xi+dx)/(2*dx))*uContorno;
	b(dim-1) -= (1.0/(dx*dx) + coefDx(xf-dx)/(2*dx))*uContorno;

	cout << "A:\n" <<  A << endl;
	cout << "b:\n" << b << "\n\n";
  	for (i = 0; i < dim-1; i++){
  		A(i+1,i+1) -= A(i,i+1)*A(i+1,i)/A(i,i);
  		b(i+1) -= b(i)*A(i+1,i)/A(i,i);
  		A(i+1,i) = 0.0;
  	}

  	b(dim-1) =  b(dim-1)/A(dim-1,dim-1);
  	for (i = dim-1; i > 0; i--){
  		b(i-1) = (b(i-1) - A(i-1,i)*b(i))/A(i-1,i-1);
  	}

  	return b;
}

VectorXd Contorno2D(double xi, double xf, double yi, double yf, double uContorno, double coefLapl,
                    funcaoReal2D coefDxDy, funcaoReal2D coefDx, funcaoReal2D coefDy, 
                    funcaoReal2D coefU, funcaoReal2D F,  int nparticoesx, int nparticoesy){
    int dim = (nparticoesy-1)*(nparticoesx-1), i, j;
    double x, y, dx = (xf - xi)/nparticoesx, dy = (yf - yi)/nparticoesy;
    VectorXd u(dim),
             b(dim);
    MatrixXd A(dim, dim);

    for(i=0; i<dim; i++){
        x = xi+dx + (i%(nparticoesx-1))*dx;
        y = yi+dy + (i/(nparticoesx-1))*dy;
        b(i) = -F(x, y);
        if(i%(nparticoesx-1) == 0){
            // b(i) -= (coefDxDy(x,y)/(4*dx*dy))*uContorno;
            b(i) -= ( coefLapl/(dx*dx) - coefDx(x,y)/(2*dx) )*uContorno;
            // b(i) -= -(coefDxDy(x,y)/(4*dx*dy))*uContorno;
            if (i/(nparticoesx-1) == 0)
                b(i) -= (coefDxDy(x,y)/(4*dx*dy))*uContorno;
        }
        if(i%(nparticoesx-1) == nparticoesx-2){
            // b(i) -= -(coefDxDy(x,y)/(4*dx*dy))*uContorno;
            b(i) -= ( coefLapl/(dx*dx) + coefDx(x,y)/(2*dx) )*uContorno;
            // b(i) -= (coefDxDy(x,y)/(4*dx*dy))*uContorno;
            if (i/(nparticoesx-1) == 0)
                b(i) -= -(coefDxDy(x,y)/(4*dx*dy))*uContorno;
        }
        if(i/(nparticoesx-1) == 0){
            // b(i) -= (coefDxDy(x,y)/(4*dx*dy))*uContorno;
            b(i) -= ( coefLapl/(dy*dy) - coefDy(x,y)/(2*dy) )*uContorno;
            // b(i) -= -(coefDxDy(x,y)/(4*dx*dy))*uContorno;
            if(i%(nparticoesx-1) == 0)
                b(i) -= -(coefDxDy(x,y)/(4*dx*dy))*uContorno;
        }
        if(i/(nparticoesx-1) == nparticoesy-2){
            // b(i) -= -(coefDxDy(x,y)/(4*dx*dy))*uContorno;
            b(i) -= ( coefLapl/(dy*dy) + coefDy(x,y)/(2*dy) )*uContorno;
            // b(i) -= +(coefDxDy(x,y)/(4*dx*dy))*uContorno;
            if(i%(nparticoesx-1) == nparticoesx-2)
                b(i) -= (coefDxDy(x,y)/(4*dx*dy))*uContorno;
        }

        A(i,i) = -2*( 1/(dx*dx) + 1/(dy*dy) );
        if(i-1 >= 0) A(i,i-1) = coefLapl/(dx*dx) - coefDx(x,y)/(2*dx);
        if(i+1 < dim) A(i,i+1) = coefLapl/(dx*dx) + coefDx(x,y)/(2*dx);
        
        if(i - nparticoesx > 0) 
            A(i,i - nparticoesx) = coefDxDy(x,y)/(4*dx*dy);
        if(i - nparticoesx+1 > 0) 
            A(i,i - nparticoesx+1) = coefLapl/(dy*dy) - coefDy(x,y)/(2*dy);
        if(i - nparticoesx+2 > 0) 
            A(i,i - nparticoesx+2) = -coefDxDy(x,y)/(4*dx*dy);

        if(i + nparticoesx-2 < dim)
            A(i,i + nparticoesx-2) = -coefDxDy(x,y)/(4*dx*dy);
        if(i + nparticoesx-1 < dim) 
            A(i,i + nparticoesx-1) = coefLapl/(dy*dy) + coefDy(x,y)/(2*dy);
        if(i + nparticoesx < dim)
            A(i,i + nparticoesx) = coefDxDy(x,y)/(4*dx*dy);

    }

    u = A.inverse() * b;

    return u;
}