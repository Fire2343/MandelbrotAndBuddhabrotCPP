#include <iostream>
#include <vector>
#include <deque>
#include <chrono>
#include <random>
#include <thread>
#include <mutex>
#include "CImg.h"

const int imgDepth = 1;
const int imgSpectrum = 3;
std::mutex vectorMutex;


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

double outsideMandelbrotOpt(double r, double i, int iterations) {
	double absSquared = r * r + i * i;
	double result = 0;
	if (absSquared * (8 * absSquared - 3) <= (double)3 / 32 - r) { //checks if inside main carinoid
		//cout << "iran" << endl;
		return result;
	}

	double zr = 0;
	double zi = 0;
	for (int iter = 0; iter < iterations; iter++) {

		double zrt = zr * zr - zi * zi + r;
		zi = 2 * zr * zi + i;
		zr = zrt;
		result = zr * zr + zi * zi;

		if (result > 2.0) {
			return result * -1;
		}
	}
	return result;
}

double F(vector<vector<double>> orbit, double rAxisMaximum, double iAxisMaximum) {
	for (int i = 0; i < orbit.size(); i++)
	{
		if ((abs(0 - orbit[i][0]) <= rAxisMaximum && abs(0 - orbit[i][1]) <= iAxisMaximum))
			return 1.0;
	}
	return 0.0;
}

double transitionProbability(double orbitSize, double mutationOrbitSize, double iterations) {
	
	double tp = (1 - (iterations - orbitSize) / iterations) / (1 - (iterations - mutationOrbitSize) / iterations);
	return tp;
}

double* mutation1(double r, double i, unsigned seed) {
	mt19937 generator(seed);
	uniform_real_distribution<double> distribution(0.0, 1.0);
	double r1 = 0.0001;
	double r2 = 0.1;
	
	double rn = r;
	double in = i;
	double phi = distribution(generator) * 3.14159265358979323846 * 2;
	r = 0.1 * exp(-log(r2 / r1) * distribution(generator));

	rn += r * cos(phi);
	in += r * sin(phi);

	double next[2] = { rn, in };
	return next;
}

vector<vector<double>> buddahBrotFunctionOpt(double r, double i, int iterations) {
	
	double absSquared = r * r + i * i;
	vector<vector<double>> result;
	
	if (absSquared * (8 * absSquared - 3) <= (double)3 / 32 - r) { //checks if inside main carinoid
		//cout << "iran" << endl;
		return result;
	}
	double zr = 0;
	double zi = 0;
	vector<double> cnumber;
	cnumber.push_back(r);
	cnumber.push_back(i);
	result.push_back(cnumber);
	for (int iter = 0; iter < iterations; iter++) {
		double zrt = zr * zr - zi * zi + r;
		zi = 2 * zr * zi + i;
		zr = zrt;
		absSquared = zr * zr + zi * zi;
		if (absSquared > 2.0) {
			break;
		}
		//if (iter > 10) {
			cnumber[0] = zr;
			cnumber[1] = zi;
			result.push_back(cnumber);
		//}
	}
	return result;
}

void mandelbrotOpt(int imageWidth, int imageHeight, double rAxisMinimum, double rAxisMaximum, double iAxisMinimum, double iAxisMaximum, int iterations) {

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
			double rn = r + rAxisMinimum;
			double in = i + iAxisMinimum;
			double result = outsideMandelbrotOpt(rn, in, iterations);
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
		cout << rPixelCounter << endl;
	}
	CImgDisplay display(mandelbrot, "", 0);
	mandelbrot.save("test.bmp");
	while (!display.is_closed()) {
		display.wait();
	}
}

void buddahbrot(double rAxisMinimum, double rAxisMaximum, double iAxisMinimum, double iAxisMaximum, int imageWidth, int imageHeight, int iterations, deque<vector<vector<double>>>* orbits, int sampleSize, int threadID) {

	
	unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
	mt19937 generator(seed);
	uniform_real_distribution<double> rDistribution(0.0, rAxisMaximum - rAxisMinimum);
	uniform_real_distribution<double> iDistribution(0.0, iAxisMaximum - iAxisMinimum);

	for (int s = 0; s < sampleSize; s++) {
		double r = rDistribution(generator);
		double i = iDistribution(generator);


		double rn = r + rAxisMinimum;
		double in = i + iAxisMinimum;
		
		
		//cout << "iran" << endl;
		//lock_guard<mutex> guard(vectorMutex);
		vectorMutex.lock();
		if ((*orbits).size() < 10000) {
			(*orbits).push_back(buddahBrotFunctionOpt(rn, in, iterations));
		}
		vectorMutex.unlock();
		

		/*if (threadID == 0) {
			cout << threadID << " " << s << endl;
		}*/
	}
}

void buddahbrotOpt(int imageWidth, int imageHeight, double rAxisMinimum, double rAxisMaximum, double iAxisMinimum, double iAxisMaximum, int iterations, int* color, CImg<double>* buddahbrot, int sampleSize) {

	unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
	mt19937 generator(seed);
	uniform_real_distribution<double> rDistribution(0.0, 4.0);
	uniform_real_distribution<double> iDistribution(0.0, 2.0);
	uniform_int_distribution<int> mDistribution(1, 5);
	uniform_int_distribution<int> alphaDistribution(0, 1);
	
	double rAbs = abs(rAxisMaximum - rAxisMinimum);
	double iAbs = abs(iAxisMaximum - iAxisMinimum);
	double rPixelValue = rAbs / imageWidth;
	double iPixelValue = iAbs / imageHeight;

	vector<int> width(imageWidth, 0);
	vector<vector<int>> matrix(imageHeight, width);
	for (int s = 0; s < sampleSize; s++) {
		double r = rDistribution(generator);
		double i = iDistribution(generator);
		double rm;
		double im;
		if (mDistribution(generator) > 1) {
			double* m = mutation1(r, i, seed);
			rm = m[0];
			im = m[1];
		}
		else {
			rm = rDistribution(generator);
			im = iDistribution(generator);
		}

		double rn = r + rAxisMinimum;
		double in = i + iAxisMinimum;

		double rnm = rm + rAxisMinimum;
		double inm = im + iAxisMinimum;

		vector<vector<double>> toColor = buddahBrotFunctionOpt(rn, in, iterations);
		vector<vector<double>> toColorM = buddahBrotFunctionOpt(rnm, inm, iterations);

		double fx = F(toColor, rAxisMaximum, rAxisMinimum);
		double fxm = F(toColorM, rAxisMaximum, rAxisMinimum);

		double tx = transitionProbability((double) toColor.size(), (double) toColorM.size(), (double) iterations);
		double txm = transitionProbability((double) toColorM.size(), (double) toColor.size(), (double) iterations);

		double alpha = (fxm * txm) / (fx * tx);
		if (alpha > alphaDistribution(generator)) {
			toColor = toColorM;
		}

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
		cout << s << endl;
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

	int selection = 0;
	while (selection == 0) {
		cout << "mandelbrot or buddahbrot? (1 for mandel, 2 for buddah)" << endl;
		cin >> selection;
	}
	int iterations = 0;
	int color[3] = {0, 0, 0};
	if (selection == 1) {
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
		cout << "insert number of iterations of the mandelbrot function for each point" << endl;
		cin >> iterations;
		mandelbrotOpt(imageWidth, imageHeight, rAxisMinimum, rAxisMaximum, iAxisMinimum, iAxisMaximum, iterations);
	}

	/*if (selection == 2) {
		cout << "insert number of iterations of the mandelbrot function for each point" << endl;
		cin >> iterations;
		int sampleSize;
		cout << "insert points sample size" << endl;
		cin >> sampleSize;
		cout << "insert r,g,b values to color buddahbrot with" << endl;
		cin >> color[0] >> color[1] >> color[2];
		CImg<double> buddahbrotImg(imageWidth, imageHeight, imgDepth, imgSpectrum, 0);
		buddahbrotOpt(imageWidth, imageHeight, rAxisMinimum, rAxisMaximum, iAxisMinimum, iAxisMaximum, iterations, color, &buddahbrotImg, sampleSize);
		
		CImgDisplay display(buddahbrotImg, "", 0);
		buddahbrotImg.save("bruddah.bmp");
		while (!display.is_closed()) {
			display.wait();
		}
	}*/
	if (selection == 2) {
		cout << "insert number of iterations of the mandelbrot function for each point" << endl;
		cin >> iterations;
		int sampleSize;
		cout << "insert points sample size" << endl;
		cin >> sampleSize;
		cout << "insert r,g,b values to color buddahbrot with" << endl;
		cin >> color[0] >> color[1] >> color[2];
		vector<thread> threads;
		int threadNumber;
		cout << "insert number of threads to use" << endl;
		cin >> threadNumber;
		CImg<double> buddahbrotImg(imageWidth, imageHeight, imgDepth, imgSpectrum, 0);

		double rAxisMaximum = 2.0;
		double rAxisMinimum = -2.0;
		double iAxisMaximum = 1.0;
		double iAxisMinimum = -1.0;
		
		
		deque<vector<vector<double>>> orbits;
		vector<int> width(imageWidth, 0);
		vector<vector<int>> matrix(imageHeight, width);
		
		double rAbs = abs(rAxisMaximum - rAxisMinimum);
		double iAbs = abs(iAxisMaximum - iAxisMinimum);
		double rPixelValue = rAbs / imageWidth;
		double iPixelValue = iAbs / imageHeight;
		

		for (int t = 0; t < threadNumber; t++) {

			thread newThread(buddahbrot, rAxisMinimum, rAxisMaximum, iAxisMinimum, iAxisMaximum, imageWidth, imageHeight, iterations, &orbits, (sampleSize / threadNumber), t);
			threads.push_back(move(newThread));
		}
		while (sampleSize > 0) {
			vectorMutex.lock();
			while (orbits.size() > 0) {
				if (orbits[0].size() < iterations + 1) {
					for (int ci = 0; ci < orbits[0].size(); ci++) {
						double r = orbits[0][ci][0];
						double i = orbits[0][ci][1];
						if (abs(0 - r) <= rAxisMaximum && abs(0 - i) <= iAxisMaximum) {
							int x = (r - rAxisMinimum) / rPixelValue;
							if (x > 0) {
								x--;
							}
							int y = (i - iAxisMinimum) / iPixelValue;
							if (y > 0) {
								y--;
							}

							//cout << x << " " << y << endl;
							matrix[y][x] = matrix[y][x] + 1;
						}
					}
				}
				orbits.pop_front();
				sampleSize--;
			}
			vectorMutex.unlock();
		}
		
		for (int t = 0; t < threadNumber; t++) {
			(threads[t]).join();
		}

		int colorInitValues[3] = {color[0], color[1], color[2]};

		for (int y = 0; y < matrix.size(); y++) {
			for (int x = 0; x < matrix[y].size(); x++) {
				color[0] = matrix[y][x] * colorInitValues[0];
				color[1] = matrix[y][x] * colorInitValues[1];
				color[2] = matrix[y][x] * colorInitValues[2];
				buddahbrotImg.draw_point(x, y, color);
			}
		}


		CImgDisplay display(buddahbrotImg, "", 0);
		buddahbrotImg.save("bruddah-naive.bmp");
		while (!display.is_closed()) {
			display.wait();
		}
	}
}
