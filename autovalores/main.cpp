#include "MN_autovalores.hpp"
#include <iostream>
using namespace std;
using namespace Eigen;

int main(int argc, char ** argv){
	MN_amb_autovalores amb(5);
	MatrixXd A(5,5);
	A << 34,3,7,7,7,
	  	3,20,3,1,10,
		7,3,64,14,14,
		7,1,14,37,10,
		7,10,14,10,53;
	cout << "Matrix entrada:\n" << A << "\n";
	VectorXd v(5);
	v << 0.052,0.482,0.884,0.679,0.357;
	
	amb.potencia(A,v,0.001);
	amb.potenciaInversa(A,v,0,0.001);
	cout << amb.autovalores() << "\n" << amb.autovetores()[4].transpose();


	return 0;
}
