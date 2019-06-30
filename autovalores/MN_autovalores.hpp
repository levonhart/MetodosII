#ifndef MN_AUTOVALORES
#define MN_AUTOVALORES
#include <Eigen/Dense>


struct MN_amb_autovalores{
  private:
	int dim;
	Eigen::VectorXd * v, w;

  public:
	MN_amb_autovalores(int n);
	Eigen::VectorXd * autovetores();
	Eigen::VectorXd autovalores();


	double potencia(Eigen::MatrixXd A, Eigen::VectorXd z, double erro);
	double potenciaInversa(Eigen::MatrixXd A, Eigen::VectorXd z, double p, double erro); // potencia( inversa de (A - p*I) )
};



#endif
