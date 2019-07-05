#ifndef MN_AUTOVALORES
#define MN_AUTOVALORES
#include <Eigen/Dense>

typedef struct parMatriz {
	Eigen::MatrixXd resultado;
	Eigen::MatrixXd acumulada;
}MP;
struct MN_amb_autovalores{
  private:
	int dim;
	Eigen::MatrixXd w;
	Eigen::MatrixXd v;

  public:
	MN_amb_autovalores(int n);
	Eigen::MatrixXd autovetores();
	Eigen::MatrixXd autovalores();


	MP potencia(Eigen::MatrixXd A, Eigen::VectorXd z, double erro);
	MP potenciaInversa(Eigen::MatrixXd A, Eigen::VectorXd z, double p, double erro); // potencia( inversa de (A - p*I) )
	MP houseHolder(Eigen::MatrixXd A);
	MP jacobi(Eigen::MatrixXd A, Eigen::MatrixXd HH, double erro);
	MP qr(Eigen::MatrixXd A, Eigen::MatrixXd Hh, double erro);


};



#endif
