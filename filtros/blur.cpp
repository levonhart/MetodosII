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
				"\t-o,\t--output <file>\t\tDetermina o local da saida. O formato "\
				"eh determinado pela extensao.\n");

}

void blur(const Mat& imagem, Mat& resultado){
	CV_Assert(imagem.depth() == CV_8U); // lança exceção se a imagem não for unsigned char
	int i, j;

	int const n = imagem.channels();

	resultado.create(imagem.size(), imagem.type());
	Mat mascara = Mat(imagem.size(), CV_16SC(n));
	for (j = 0; j < imagem.rows; j++) {
		for (i = 0; i < n*imagem.cols; i++) {
			if (i<n) {
				//foward x
				mascara.ptr<int16_t>(j)[i] = (imagem.ptr<uchar>(j)[i+2*n] - 2*imagem.ptr<uchar>(j)[i+n] + imagem.ptr<uchar>(j)[i])/4;
			} else if (i>=n*(imagem.cols-1)) {
				//backward x
				mascara.ptr<int16_t>(j)[i] = (imagem.ptr<uchar>(j)[i] - 2*imagem.ptr<uchar>(j)[i-n] + imagem.ptr<uchar>(j)[i-2*n])/4;
			} else{
				//central x
				mascara.ptr<int16_t>(j)[i] = (imagem.ptr<uchar>(j)[i+n] - 2*imagem.ptr<uchar>(j)[i] + imagem.ptr<uchar>(j)[i-n])/4;
			}
			if (j==0) {
				//foward y
				mascara.ptr<int16_t>(j)[i] += (imagem.ptr<uchar>(j+2)[i] - 2*imagem.ptr<uchar>(j+1)[i] + imagem.ptr<uchar>(j)[i])/4;
			} else if (j==imagem.rows-1) {
				//backward y
				mascara.ptr<int16_t>(j)[i] += (imagem.ptr<uchar>(j)[i] - 2*imagem.ptr<uchar>(j-1)[i-n] + imagem.ptr<uchar>(j-2)[i])/4;
			} else{
				//central y
				mascara.ptr<int16_t>(j)[i] += (imagem.ptr<uchar>(j+1)[i] - 2*imagem.ptr<uchar>(j)[i] + imagem.ptr<uchar>(j-1)[i])/4;
			}
			resultado.ptr<uchar>(j)[i] = saturate_cast<uchar>(imagem.ptr<uchar>(j)[i] + mascara.ptr<int16_t>(j)[i]);
			//mascara.ptr<int16_t>(j)[i] = saturate_cast<uchar>(mascara.ptr<int16_t>(j)[i]+sizeof(uchar)/2);

		}	
	}
	mascara.convertTo(mascara,CV_8UC(n),2,0);
	namedWindow("Máscara",WINDOW_AUTOSIZE);
	imshow("Máscara",mascara);

}

int main(int argc, char *argv[]){
	String caminho;
	String saida = "saida.bmp";
	int i;
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
			}
		}
	} else {
		help();
		return -2;
	}
    Mat	imagem, blur_imagem;
	
	imagem = imread(caminho,IMREAD_UNCHANGED);
	if(imagem.empty()){
		printf("Nao foi possivel ler a imagem.\n");
		return -1;
	}
	blur(imagem, blur_imagem);

	namedWindow(caminho,WINDOW_AUTOSIZE);
	imshow(caminho,imagem);
	namedWindow("saida",WINDOW_AUTOSIZE);
	imshow("saida",blur_imagem);
	imwrite(saida,blur_imagem);
	waitKey(0);
	return 0;
}
