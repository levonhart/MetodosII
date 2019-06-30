#include "MN_autovalores.hpp"
#include <cmath>
#include <iostream>

using namespace std;
using namespace Eigen;

MN_amb_autovalores::MN_amb_autovalores(int n){
	this->dim = n;
	this->w.resize(n,1);
	this->v = (VectorXd*)malloc(n*sizeof(VectorXd));
}

VectorXd * MN_amb_autovalores::autovetores(){
	return v;
}
VectorXd MN_amb_autovalores::autovalores(){
	return w;
}


double MN_amb_autovalores::potencia(MatrixXd A, VectorXd z, double erro){
	VectorXd z_ant(dim), z_act(dim);
	double w_ant, w_apx, er_rel = +INFINITY;
	z_act= (1.0/z.lpNorm<Infinity>())*z;
	w_apx = z_act.dot(A*z_act) / z_act.dot(z_act);
	cout << "w = " << w_apx << "	z = " << z_act.transpose() << "	epsilon = "<< er_rel << "\n";
	while(er_rel > erro){
		w_ant = w_apx;
		z_ant = z_act;
		z_act = A*z_ant;
		z_act = (1.0/z_act.lpNorm<Infinity>())*z_act;
		w_apx = z_act.dot(A*z_act) / z_act.dot(z_act);
		er_rel = abs( (w_apx - w_ant)/w_apx );
		cout << "w = " << w_apx << "	z = " << z_act.transpose() << "	epsilon = "<< er_rel << "\n";
	}

	this->w(dim-1) = w_apx;
	this->v[dim-1] = z_act;
	
	return w_apx;
}

double MN_amb_autovalores::potenciaInversa(MatrixXd A, VectorXd z, double p, double erro){
	double tmp1 = w(dim-1);
	VectorXd tmp2(dim);
	tmp2 = v[dim-1];


	MatrixXd B = (A - p*MatrixXd::Identity(dim,dim));
	B = B.inverse();
	potencia(B,z,erro);

	w(0) = 1/w(dim-1)+p;
	w(dim-1) = tmp1;
	v[0] = v[dim-1];
	v[dim-1] = tmp2;

	return w(0);
}
