#include <iostream>
#include "CImg.h"

using namespace std;
using namespace cimg_library;

double* squareOfComplex(double number[2]) {
	double r = number[0];
	double i = number[1];
	double result[2];
	result[0] = (r * r) - (i * i);
	result[1] = 2 * (r * i);
	return result;
}

double complexAbs(double number[2]) {
	double r = number[0];
	double i = number[1];
	double result = sqrt((r * r) + (i * i));
	return result;
}

double outsideMandelbrot(double number[2], int iterations) {
	double z[2] = {0.0, 0.0};
	double result;
	for (int iter = 0; iter < iterations; iter++) {
		double* square = squareOfComplex(z);
		z[0] = square[0] + number[0];
		z[1] = square[1] + number[1];
		result = complexAbs(z);
		if (result > 2.0) {
			return result * -1;
		}
	}
	return result;
}

int main(){
	int imageWidth;
	int imageHeight;
	cout << "please input the width and height of the desired image, in pixels" << endl;
	cin >> imageWidth >> imageHeight;
	double rAxisMinimum;
	double rAxisMaximum;
	double iAxisMinimum;
	double iAxisMaximum;
	cout << "please input the minimum and maximum real axis coordinates of the desired image" << endl;
	cin >> rAxisMinimum >> rAxisMaximum;
	double rAbs = abs(rAxisMaximum - rAxisMinimum);
	cout << "please input the minimum and maximum imaginary axis coordinates of the desired image" << endl;
	cin >> iAxisMinimum >> iAxisMaximum;
	double iAbs = abs(iAxisMaximum - iAxisMinimum);
	CImg<double> mandelbrot(imageWidth, imageHeight, 1, 3, 0);
	double color[3];
	color[0] = 0;
	color[1] = 0;
	color[2] = 0;
	double rPixelValue = rAbs / imageWidth;
	double iPixelValue = iAbs / imageHeight;
	int rPixelCounter = 0;
	int iPixelCounter = 0;
	for (double r = 0; r < rAbs; r += rPixelValue) {
		for (double i = 0; i < iAbs; i += iPixelValue) {
			double number[2] = {r + rAxisMinimum, i + iAxisMinimum};
			double result = outsideMandelbrot(number, 100);
			if (result < 0) {
				color[0] = 40 / abs(result);
				color[2] = 50 / abs(result);
				mandelbrot.draw_point(rPixelCounter, imageHeight - iPixelCounter, color);				
			}
			iPixelCounter++;
		}
		iPixelCounter = 0;
		rPixelCounter++;
	}
	CImgDisplay display(mandelbrot, "");
	while (!display.is_closed()) {
		display.wait();
	}

}
