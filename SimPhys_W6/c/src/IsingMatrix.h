/*
 * IsingMatrix.h
 *
 *  Created on: 01.02.2013
 */

#ifndef ISINGMATRIX_H_
#define ISINGMATRIX_H_

#include <iostream>
#include <gsl/gsl_rng.h> //libgsl0-dev

using namespace std;

class IsingMatrix {
public:
	IsingMatrix(int n, int m, unsigned long int seed = 42, double J = 1, double H = 0);
	~IsingMatrix();
	void init();
	void flip();
	void unflip();
	double getE();
	double getM();
	int& operator() (int x, int y);
	int* operator[] (int i);

private:
	int PlusX(int i);
	int PlusY(int i);
	int MinusX(int i);
	int MinusY(int i);
	int dx,dy,s;
	double J, H, E, M;
	int *matrix;
	gsl_rng *rnd;
	int iold;
	double Eold, Mold;
};

#endif /* ISINGMATRIX_H_ */
