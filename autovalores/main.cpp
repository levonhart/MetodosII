#include "MN_autovalores.hpp"
#include <iostream>
#include <iomanip>
using namespace std;
using namespace Eigen;

#define EPSILON 0.01e-6

int main(int argc, char ** argv){
	MN_amb_autovalores amb(5);
	MP resul;
	MatrixXd A(5,5), U(5,5), S(5,5), V(5,5), B;
	A << 34,3,7,7,7,
	  	3,20,3,1,10,
		7,3,64,14,14,
		7,1,14,37,10,
		7,10,14,10,53;

	cout << "=====================================================\n";
	cout << "==================== QUESTAO 01 ==================\n";
	cout << "=====================================================\n";

	cout << "Matrix entrada:\n\n" << A << "\n\n";
	VectorXd v(5);
	v << 0.052,0.482,0.884,0.679,0.357;
	cout << "vetor inicial metodo da potencia (transposto):\n\n";
	cout << v.transpose() << "\n\n";
	
	cout<<"metodo da potencia\n\n";
	resul = amb.potencia(A,v,EPSILON);
	cout << "autovalor:\t" << resul.resultado << "\n\n";
	cout << "autovetor:\t" << resul.acumulada.transpose()<<"\n\n";

	cout<<"metodo da potencia inversa\n\n";
	resul = amb.potenciaInversa(A,v,0,EPSILON);
	cout << "autovalor:\t" << resul.resultado << "\n\n";
	cout << "autovetor:\t" << resul.acumulada.transpose()<<"\n\n";

	cout << "=====================================================\n\n";
	




	cout << "=====================================================\n";
	cout << "==================== QUESTAO 02 ==================\n";
	cout << "=====================================================\n";
	
	cout<<"metodo de jacobi\n\n";
	amb.jacobi(A*A.transpose(),MatrixXd::Identity(5,5),EPSILON);
	U = amb.autovetores();
	S.block = amb.autovalores().cwiseAbs().cwiseSqrt();

	cout << "autovalores AAt:\n" << amb.autovalores()<<"\n\n";
	cout << "Sigma:\n" << S << "\n\n";
	cout << "U:\n" << U <<"\n\n";

	amb.jacobi(A.transpose()*A,MatrixXd::Identity(5,5),EPSILON);
	V = amb.autovetores(); 
	S = amb.autovalores().cwiseAbs().cwiseSqrt();
	cout << "autovalores AtA:\n" << amb.autovalores()<<"\n\n";
	cout << "Sigma:\n" << S << "\n\n";
	cout << "V:\n" << V <<"\n\n";


	cout << "=====================================================\n";
	cout << "USVt:\n" << U*S*V.transpose() <<"\n\n";
 	cout << "=====================================================\n";


	cout<<"metodo QR\n\n";
	amb.qr(A*A.transpose(),MatrixXd::Identity(5,5),EPSILON);
	U = amb.autovetores();
	S = amb.autovalores().cwiseAbs().cwiseSqrt();
	cout << "autovalores AAt:\n" << amb.autovalores()<<"\n\n";
	cout << "Sigma:\n" << S << "\n\n";
	cout << "U:\n" << U <<"\n\n";

	amb.qr(A.transpose()*A,MatrixXd::Identity(5,5),EPSILON);
	V = amb.autovetores(); 
	S = amb.autovalores().cwiseAbs().cwiseSqrt();
	cout << "autovalores AtA:\n" << amb.autovalores()<<"\n\n";
	cout << "Sigma:\n" << S << "\n\n";
	cout << "V:\n" << V <<"\n\n";

	cout << "=====================================================\n";
	cout << "USVt:\n" << U*S*V.transpose() <<"\n\n";


	// cout<<"matriz de house-holder\n\n";
	// resul = amb.houseHolder(A);
	// cout << "matriz Hh:\n" << resul.resultado << "\n\n";
	// cout << "acumulada:\n" << resul.acumulada <<"\n\n";

	// cout<<"metodo de jacobi com house-holder\n\n";
	// amb.jacobi(resul.resultado,resul.acumulada,EPSILON);
	// cout << "autovalores:\n" << amb.autovalores() << "\n\n";
	// cout << "autovetores:\n" << amb.autovetores() <<"\n\n";


	// cout<<"metodo QR com house-holder\n\n";
	// amb.qr(resul.resultado,resul.acumulada,EPSILON);
	// cout << "autovalores:\n" << amb.autovalores() << "\n\n";
	// cout << "autovetores:\n" << amb.autovetores() <<"\n\n";

	cout << "=====================================================\n\n";

	



	cout << "=====================================================\n";
	cout << "==================== QUESTAO 03 ==================\n";
	cout << "=====================================================\n";
	
	A = MatrixXd::Zero(3,5);
	U = MatrixXd::Zero(3,3);
	V = MatrixXd::Zero(5,5);
	S = MatrixXd::Zero(3,5);

	A << 34,3,7,7,7,
	  	3,20,3,1,10,
		7,3,64,14,14;

	cout << "Matrix entrada:\n\n" << A << "\n\n";
	
	cout<<"metodo de jacobi\n\n";
	B = A*A.transpose();
	cout << "AAt:\n" << B<<"\n\n";
	amb.jacobi(B,MatrixXd::Identity(3,3),EPSILON);


	U = amb.autovetores();
	S.block<3,3>(0,0) = amb.autovalores().cwiseAbs().cwiseSqrt().block<3,3>(0,0);
	cout << "autovalores AAt:\n" << amb.autovalores()<<"\n\n";
	cout << "Sigma:\n" << S << "\n\n";
	cout << "U:\n" << U <<"\n\n";

	B = A.transpose()*A;
	cout << "AtA:\n" << B<<"\n\n";
	amb.jacobi(B,MatrixXd::Identity(5,5),EPSILON);
	V = amb.autovetores();
	S.block<3,3>(0,0) = amb.autovalores().cwiseAbs().cwiseSqrt().block<3,3>(0,0);
	cout << "autovalores AtA:\n" << amb.autovalores()<<"\n\n";
	cout << "Sigma:\n" << S << "\n\n";
	cout << "V:\n" << V <<"\n\n";


	cout << "=====================================================\n";
	cout << "USVt:\n" << U*S*V.transpose() <<"\n\n";
 	cout << "=====================================================\n";


	cout<<"metodo QR\n\n";

	B = A*A.transpose();
	cout << "AAt:\n" << B<<"\n\n";
	amb.qr(B,MatrixXd::Identity(3,3),EPSILON);
	U = amb.autovetores();
	S.block<3,3>(0,0) = amb.autovalores().cwiseAbs().cwiseSqrt().block<3,3>(0,0);
	cout << "autovalores AAt:\n" << amb.autovalores()<<"\n\n";
	cout << "Sigma:\n" << S << "\n\n";
	cout << "U:\n" << U <<"\n\n";

	B = A.transpose()*A;
	cout << "AtA:\n" << B<<"\n\n";
	amb.qr(B,MatrixXd::Identity(5,5),EPSILON);
	V = amb.autovetores();
	S.block<3,3>(0,0) = amb.autovalores().cwiseAbs().cwiseSqrt().block<3,3>(0,0);
	cout << "autovalores AtA:\n" << amb.autovalores()<<"\n\n";
	cout << "Sigma:\n" << S << "\n\n";
	cout << "V:\n" << V <<"\n\n";

	cout << "=====================================================\n";
	cout << "USVt:\n" << U*S*V.transpose() <<"\n\n";


	return 0;
}
