#include <iostream>
#include <vector>
#include <random>
#include <thread>
#include "CImg.h"

const int imgDepth = 1;
const int imgSpectrum = 3;


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

vector<vector<double>> buddahBrotFunction(double number[2], int iterations) {
	double z[2] = {0.0 , 0.0};
	double abs;
	vector<vector<double>> result;
	vector<double> cnumber;
	cnumber.push_back(number[0]);
	cnumber.push_back(number[1]);
	result.push_back(cnumber);
	for (int iter = 0; iter < iterations; iter++) {
		double* square = squareOfComplex(z);
		z[0] = square[0] + number[0];
		z[1] = square[1] + number[1];
		abs = complexAbs(z);
		if (abs > 2.0) {
			break;
		}
		//cout << "z is" << z[0] << " " << z[1] << endl;
		cnumber[0] = z[0];
		cnumber[1] = z[1];
 		result.push_back(cnumber);
		//cout << "result is" << result[result.size() - 1][0] << " " << result[result.size() - 1][1] << endl;
	}
	return result;
}

void mandelbrot(int imageWidth, int imageHeight, double rAxisMinimum, double rAxisMaximum, double iAxisMinimum, double iAxisMaximum, int iterations) {
	
	double rAbs = abs(rAxisMaximum - rAxisMinimum);
	double iAbs = abs(iAxisMaximum - iAxisMinimum);
	CImg<double> mandelbrot(imageWidth, imageHeight, imgDepth, imgSpectrum, 0);
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
			double number[2] = { r + rAxisMinimum, i + iAxisMinimum };
			double result = outsideMandelbrot(number, iterations);
			if (result < 0) {
				color[0] = 255 / abs(result);
				color[1] = 0 / abs(result);
				color[2] = 0 / abs(result);
				mandelbrot.draw_point(rPixelCounter, imageHeight - iPixelCounter, color);
			}
			iPixelCounter++;
		}
		iPixelCounter = 0;
		rPixelCounter++;
	}
	CImgDisplay display(mandelbrot, "", 0);
	mandelbrot.save("test.bmp");
	while (!display.is_closed()) {
		display.wait();
	}
}

void buddahbrot(int imageWidth, int imageHeight, double rAxisMinimum, double rAxisMaximum, double iAxisMinimum, double iAxisMaximum, int iterations, int* color, CImg<double>* buddahbrot) {
	
	double rAbs = abs(rAxisMaximum - rAxisMinimum);
	double iAbs = abs(iAxisMaximum - iAxisMinimum);
	double rPixelValue = rAbs / imageWidth;
	double iPixelValue = iAbs / imageHeight;
	
	vector<int> width(imageWidth, 0);
	vector<vector<int>> matrix(imageHeight, width);
	for (int yPixel = 0; yPixel < imageHeight; yPixel++) {
		for (int xPixel = 0; xPixel < imageWidth; xPixel++) {
			double r = xPixel * rPixelValue;
			double i = yPixel * iPixelValue;
			double number[2] = { r + rAxisMinimum, i + iAxisMinimum };
			vector<vector<double>> toColor = buddahBrotFunction(number, iterations);
			if (toColor.size() < iterations + 1) {
				for (int i = 0; i < toColor.size(); i++) {
					if (abs(0 - toColor[i][0]) <= rAxisMaximum && abs(0 - toColor[i][1]) <= iAxisMaximum) {
						int x = (toColor[i][0] - rAxisMinimum) / rPixelValue;
						if (x > 0) {
							x--;
						}
						int y = (toColor[i][1] - iAxisMinimum) / iPixelValue;
						if (y > 0) {
							y--;
						}
						//cout << "i is" << i << endl;
						//cout << toColor[i][0] << " " << toColor[i][1] << endl;
						//cout << x << " " << y << endl;
						matrix[y][x] = matrix[y][x] + 1;
					}
				}
			}
		}
		//cout << yPixel << endl;
	}
	
	vector <int> colorInitValues;
	colorInitValues.push_back(color[0]);
	colorInitValues.push_back(color[1]);
	colorInitValues.push_back(color[2]);

	for (int y = 0; y < matrix.size(); y++) {
		for (int x = 0; x < matrix[y].size(); x++) {
			color[0] = matrix[y][x] * colorInitValues[0];
			color[1] = matrix[y][x] * colorInitValues[1];
			color[2] = matrix[y][x] * colorInitValues[2];
			(*buddahbrot).draw_point(x, y, color);
		}
	}
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

	int selection = 0;
	while (selection == 0) {
		cout << "mandelbrot or buddahbrot? (1 for mandel, 2 for buddah, 3 for buddahRGB)" << endl;
		cin >> selection;
	}
	int iterations = 0;
	int color[] = { 0, 0, 0 };
	if (selection == 1) {
		cout << "insert number of iterations of the mandelbrot function for each point" << endl;
		cin >> iterations;
		mandelbrot(imageWidth, imageHeight, rAxisMinimum, rAxisMaximum, iAxisMinimum, iAxisMaximum, iterations);
	}

	if (selection == 2) {
		cout << "insert number of iterations of the mandelbrot function for each point" << endl;
		cin >> iterations;
		cout << "insert r,g,b values to color buddahbrot with" << endl;
		cin >> color[0] >> color[1] >> color[2];
		CImg<double> buddahbrotImg(imageWidth, imageHeight, imgDepth, imgSpectrum, 0);
		buddahbrot(imageWidth, imageHeight, rAxisMinimum, rAxisMaximum, iAxisMinimum, iAxisMaximum, iterations, color, &buddahbrotImg);
		
		CImgDisplay display(buddahbrotImg, "", 0);
		buddahbrotImg.save("bruddah.bmp");
		while (!display.is_closed()) {
			display.wait();
		}
	}
	/*if (selection == 3) {
		iterations = 100;
		color[0] = 255;
		color[1] = 0;
		color[2] = 0;
		buddahbrot(imageWidth, imageHeight, rAxisMinimum, rAxisMaximum, iAxisMinimum, iAxisMaximum, iterations, color);
	}*/
}
