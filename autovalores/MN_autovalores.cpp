#include "MN_autovalores.hpp"
#include <cmath>
#include <iostream>

using namespace std;
using namespace Eigen;



MN_amb_autovalores::MN_amb_autovalores(int n){
	this->dim = n;
	this->w.Zero(n,n);
	this->v = MatrixXd(n, n);
}

MatrixXd MN_amb_autovalores::autovetores(){
	return v;
}
MatrixXd MN_amb_autovalores::autovalores(){
	return w;
}


MP MN_amb_autovalores::potencia(MatrixXd A, VectorXd z, double erro){
	VectorXd z_ant(dim), z_act(dim);
	double w_ant, w_apx, er_rel = +INFINITY;
	z_act= (1.0/z.lpNorm<Infinity>())*z;
	w_apx = z_act.dot(A*z_act) / z_act.dot(z_act);
	// cout << "w = " << w_apx << "	z = " << z_act.transpose() << "	epsilon = "<< er_rel << "\n";
	while(er_rel > erro){
		w_ant = w_apx;
		z_ant = z_act;
		z_act = A*z_ant;
		z_act = (1.0/z_act.lpNorm<Infinity>())*z_act;
		w_apx = z_act.dot(A*z_act) / z_act.dot(z_act);
		er_rel = abs( (w_apx - w_ant)/w_apx );
		// cout << "w = " << w_apx << "	z = " << z_act.transpose() << "	epsilon = "<< er_rel << "\n";
	}


	VectorXd lambda(1);
	lambda(0) = w_apx;
	return { lambda, z_act};
}

MP MN_amb_autovalores::potenciaInversa(MatrixXd A, VectorXd z, double p, double erro){
	MP resul;

	MatrixXd B = (A - p*MatrixXd::Identity(dim,dim));
	B = B.inverse();
	resul = potencia(B,z,erro);

	resul.resultado(0) = 1.0/resul.resultado(0) + p;
	return {resul};
}

MP MN_amb_autovalores::houseHolder(MatrixXd A){
	int n = A.rows();



	double alfa, r;
	MatrixXd Hac = MatrixXd::Identity(n, n), 
			 Hh;
	VectorXd v = VectorXd::Zero(n), vl = VectorXd::Zero(n),
			 N(n);
	for (int i = 0; i < n-2; i++){
		

		for (int k = i+1; k < n; k++)
			v(k) = A(k,i);

		alfa = v.norm();

		if(alfa > 0.1e-8){
			if(v(i+1) > 0)
				alfa = -alfa;

			r = sqrt(0.5*(alfa*alfa - v(i+1)*alfa ));

			vl(i+1) = alfa;
			N = v - vl;

			N = (0.5/r)*N;
			Hh = MatrixXd::Identity(n, n) - 2*N*N.transpose();



			A = Hh * A *Hh;
			Hac = Hac * Hh;
		}
			v.setZero(); vl.setZero(); N.setZero();
	}

	return {A, Hh};
}

double error(Ref<MatrixXd> A){
	double soma = 0;
	for (int i = 0; i < A.rows(); ++i)
		for (int j = 0; j < A.cols(); ++j)
			if(i!=j) soma += pow(A(i,j),2);
	return sqrt(soma);
}


MP MN_amb_autovalores::jacobi(MatrixXd A, MatrixXd HH, double erro){
	int n = A.rows();
	double theta;

	MatrixXd Pij = MatrixXd::Identity(n,n), J(n, n), P(n, n);

	J = HH;


	int i = 0;
	while(error(A)> erro && i < 1000000){
		P = MatrixXd::Identity(n, n);

		for (int j = 0; j < n-1; j++){
			for (int k = j+1; k < n; k++){

				theta = (fabs(A(j, j) - A(k, k))  < erro) ? (M_PI/4) : ( (0.5)*atan( (2*A(k,j)/ (A(j,j) - A(k,k)) ) ) ); 

				Pij(k,k) = cos(theta);
				Pij(j,j) = cos(theta);
				Pij(k,j) = sin(theta);
				Pij(j,k) = -sin(theta);




				A= Pij.transpose()*A*Pij;
				P *= Pij;
				
				Pij = MatrixXd::Identity(n, n);
			}
		}
		J *= P;
		i++;
	}	
	this->w = A;
	this->v = J;
	return {A,J};
}


MP MN_amb_autovalores::qr(Eigen::MatrixXd A, Eigen::MatrixXd Hh, double erro){
	int n = A.rows();
	double theta;

	MatrixXd Q = MatrixXd::Identity(n, n), R = A, P(n, n),
			Pt;

	P = Hh;

	int i=0;
	while(error(A) > erro && i < 1000000){

		for (int j = 0; j < n-1; j++){
			for (int k = (j+1); k < n; k++){

				Pt = MatrixXd::Identity(n, n);

				theta = (fabs(A(j, j)) < erro) ? (M_PI/4) : ( atan( ( A(k,j)/A(j,j) ) ) ); 

				Pt(k,k) = cos(theta);
				Pt(j,j) = cos(theta);
				Pt(k,j) = -sin(theta);
				Pt(j,k) = sin(theta);


				R  = Pt*R;
				Q = Pt*Q;			
			}
		}
		

		Q.transposeInPlace();
		
		A = R*Q;
		P *= Q;

		Q = MatrixXd::Identity(n, n);
		R = A;
		i++;
	}

	this->w = A;
	this->v = P;
	return {A,P};
}
