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
std::condition_variable cv;
bool ready = false;
bool processed = false;


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

int outsideMandelbrotOpt(double r, double i, int iterations) {
	double absSquared = r * r + i * i;
	double result = 0;
	int iter = 0;
	if (absSquared * (8 * absSquared - 3) <= (double)3 / 32 - r) { //checks if inside main carinoid
		return result;
	}

	double zr = 0;
	double zi = 0;
	while (iter < iterations) {

		double zrt = zr * zr - zi * zi + r;
		zi = 2 * zr * zi + i;
		zr = zrt;
		result = sqrt(zr * zr + zi * zi);
		iter++;

		if (result > 2.0) {
			break;
		}
	}
	return iter;
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
		cnumber[0] = zr;
		cnumber[1] = zi;
		result.push_back(cnumber);
	}
	return result;
}

void mandelbrotOpt(int imageWidth, int imageHeight, double rAxisCenter, double iAxisCenter, int iterations, double zoomFactor) {

	
	double imgWidthToHeightRatio = ((double)imageWidth) / ((double)imageHeight);
	double rHalfWidth = (1.0  * imgWidthToHeightRatio) / zoomFactor;
	double iHalfWidth = 1.0 / zoomFactor;

	double rAxisMinimum = rAxisCenter - rHalfWidth;
	double rAxisMaximum = rAxisCenter + rHalfWidth;
	double iAxisMinimum = iAxisCenter - iHalfWidth;
	double iAxisMaximum = iAxisCenter + iHalfWidth;


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
			int iters = outsideMandelbrotOpt(rn, in, iterations);
			if (iters < iterations - 1) {
				color[0] = 0.3 * iters;
				color[1] = 0.5 * iters;
				color[2] = 1.0 * iters;
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
		
		
		if ((*orbits).size() < 10000) {
			(*orbits).push_back(buddahBrotFunctionOpt(rn, in, iterations));
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
		double rAxisCenter;
		double iAxisCenter;
		double zoomFactor;
		cout << "please input the center real axis coordinates of the desired image" << endl;
		cin >> rAxisCenter;
		cout << "please input the center imaginary axis coordinates of the desired image" << endl;
		cin >> iAxisCenter;
		cout << "insert zoom factor (>=1)" << endl;
		cin >> zoomFactor;
		cout << "insert number of iterations of the mandelbrot function for each point" << endl;
		cin >> iterations;
		mandelbrotOpt(imageWidth, imageHeight, rAxisCenter, iAxisCenter, iterations, zoomFactor);
	}

	if (selection == 2) {
		cout << "insert number of iterations of the mandelbrot function for each point" << endl;
		cin >> iterations;
		int sampleSize;
		cout << "insert points sample size" << endl;
		cin >> sampleSize;
		cout << "insert r,g,b values to color buddahbrot with" << endl;
		cin >> color[0] >> color[1] >> color[2];
		vector<thread> threads;

		auto start = std::chrono::system_clock::now();
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
		


		thread newThread(buddahbrot, rAxisMinimum, rAxisMaximum, iAxisMinimum, iAxisMaximum, imageWidth, imageHeight, iterations, &orbits, sampleSize, 0);
		threads.push_back(move(newThread));
		
		int sampleCounter = 0;
		int initSampleSize = sampleSize;
		while (sampleSize > 0) {
			if (sampleCounter >= initSampleSize - 10000) {
				if (orbits.size() == 0) {
					if (orbits.size() == 0) {
						break;
					}
				}
			}
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

							matrix[y][x] = matrix[y][x] + 1;
						}
					}
				}
				orbits.pop_front();
				sampleSize--;
				sampleCounter++;
				if ((sampleCounter >= 100000 && sampleCounter % 100000 == 0) || sampleCounter >= initSampleSize - 10000) {
					cout << "sample counter is: " << sampleCounter << endl;
					cout << "sample size is: " << sampleSize << endl;
					cout << "orbits size is: " << orbits.size() << endl;
				}
			}
		}
		
		(threads[0]).join();
		
		int colorInitValues[3] = {color[0], color[1], color[2]};

		for (int y = 0; y < matrix.size(); y++) {
			for (int x = 0; x < matrix[y].size(); x++) {
				color[0] = matrix[y][x] * colorInitValues[0];
				color[1] = matrix[y][x] * colorInitValues[1];
				color[2] = matrix[y][x] * colorInitValues[2];
				buddahbrotImg.draw_point(x, y, color);
			}
		}

		auto end = std::chrono::system_clock::now();
		std::chrono::duration<double> diff = end - start;
		cout << "program took " << diff.count() << " seconds to run" << endl;
		
		CImgDisplay display(buddahbrotImg, "", 0);
		buddahbrotImg.save("buddah-naive.bmp");
		while (!display.is_closed()) {
			display.wait();
		}
	}
}
