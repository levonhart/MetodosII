#include <iostream>
#include <opencv2/core.hpp>
#include <opencv2/imgcodecs.hpp>
#include <opencv2/highgui.hpp>
#include <string>

using namespace std;
using namespace cv;

static void help(void){
		printf("\nUso: contorno [OPTIONS] PATH\n\n"\
				"\'contorno\' realca arestas e limites de objetos de uma imagem em "\
				"PATH e retorna uma imagem em B&W com os contornos\n\n"\
				"Exemplo:\ncontorno imagem/exemplo.png\n"\
				"Opcoes:\n"\
				"\t-h,\t--help       \t\tMostra essa informacao.\n"\
				"\t-l,\t--limite <arg>\t\tDetermina o limite de sensibilidade.\n"\
				"\t-o,\t--output <file>\t\tDetermina o local da saida. O formato "\
				"eh determinado pela extensao.\n");

}

void realcaContornos(const Mat& imagem, Mat& resultado, uchar limite){
	CV_Assert(imagem.depth() == CV_8U); // lança exceção se a imagem não for unsigned char
	int i, j;

	resultado.create(imagem.size(), imagem.type());
	for (j = 0; j < imagem.rows; j++) {
		for (i = 0; i < imagem.cols; i++) {
			if (i==0) {
				//foward
				resultado.ptr<uchar>(j)[i] = saturate_cast<uchar>(2*abs((imagem.ptr<uchar>(j)[i+1] - imagem.ptr<uchar>(j)[i])));
			} else if (i==imagem.cols-1) {
				//backward
				resultado.ptr<uchar>(j)[i] = saturate_cast<uchar>(abs(2*(imagem.ptr<uchar>(j)[i] - imagem.ptr<uchar>(j)[i-1])));
			} else{
				//central
				resultado.ptr<uchar>(j)[i] = saturate_cast<uchar>(abs(imagem.ptr<uchar>(j)[i+1] - imagem.ptr<uchar>(j)[i-1]));
			}
			if (limite > 0 && limite < 255) {
				if(resultado.ptr<uchar>(j)[i] > limite)
					resultado.ptr<uchar>(j)[i] = 255;
				else resultado.ptr<uchar>(j)[i] = 0;
			}
		}	
	}

}

int main(int argc, char *argv[]){
	String caminho;
	String saida = "saida.bmp";
	int i;
		uchar limite = 0;
	if (argc > 1) {
		for (i = 1; i < argc; i++) {
			if (argv[i][0] != '-' && caminho.empty()) {
				caminho = argv[i];
			} else if (!strcmp(argv[i], "--help") || !strcmp(argv[i], "-h")) {
				help();
				return 0;
			} else if (!strcmp(argv[i], "-o") || !strcmp(argv[i], "--output")) {
				i++;
					saida = argv[i];
			} else if (!strcmp(argv[i], "-l") || !strcmp(argv[i], "--limite")) {
				i++;
				limite = atoi(argv[i]);
			}
		}
		// if (!strcmp(argv[1],"--help") || !strcmp(argv[1],"-h")) {
		//     help() ;
		//     return 0;
		// }
		// caminho = argv[1];
	} else {
		help();
		return -2;
	}
    Mat	imagem, contornos;
	
	imagem = imread(caminho,IMREAD_GRAYSCALE);
	if(imagem.empty()){
		printf("Nao foi possivel ler a imagem.\n");
		return -1;
	}
	realcaContornos(imagem, contornos, limite);

	namedWindow(caminho,WINDOW_AUTOSIZE);
	imshow(caminho,imagem);
	namedWindow("saida",WINDOW_AUTOSIZE);
	imshow("saida",contornos);
	imwrite(saida,contornos);
	waitKey(0);
	return 0;
}
