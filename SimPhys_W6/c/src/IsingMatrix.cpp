/*
 * IsingMatrix.cpp
 *
 *  Created on: 01.02.2013
 */

#include "IsingMatrix.h"

IsingMatrix::IsingMatrix(int n, int m, unsigned long int seed, double J, double H) : dx(n), dy(m), s(n*m), J(J), H(H), matrix(new int [m*n]){
	rnd = gsl_rng_alloc(gsl_rng_taus);
	gsl_rng_set(rnd, seed);
	for (int i = 0; i < s; ++i) {
		matrix[i] = 0;
	}
	iold = 0;
	Eold = 0;
	Mold = 0;
}

IsingMatrix::~IsingMatrix() {
	delete[] matrix;
	gsl_rng_free(rnd);
}

void IsingMatrix::init() {
	for (int i = 0; i < s; ++i) {
		matrix[i] = 2*gsl_rng_uniform_int(rnd,2)-1;
	}

	E = 0, M = 0;
	for (int i = 0; i < s; ++i) {
		//cout << i << " " << iPx << " " << iPy << endl;
		E -= matrix[i]*(J*(matrix[PlusX(i)]+matrix[PlusY(i)])+H);
		M += matrix[i];
	}
}

int IsingMatrix::PlusX(int i) {
	return ((i+1) % dx == 0)? i+1-dx : i+1;
}

int IsingMatrix::PlusY(int i) {
	return (i+dx > s-1)	 	? i+dx-s : i+dx;
}

int IsingMatrix::MinusX(int i) {
	return (i % dx == 0) 	? i-1+dx : i-1;
}

int IsingMatrix::MinusY(int i) {
	return (i-dx < 0)	 	? i-dx+s : i-dx;
}

void IsingMatrix::flip() {
	int i = gsl_rng_uniform_int(rnd,s);

	iold = i;
	Eold = E;
	Mold = M;

	matrix[i] *= -1;
	E -= 2*matrix[i]*(J*(matrix[PlusX(i)]+matrix[MinusX(i)]+matrix[PlusY(i)]+matrix[MinusY(i)])+H);
	M += 2*matrix[i];
}

void IsingMatrix::unflip() {
	matrix[iold] *= -1;
	E = Eold;
	M = Mold;
}

double IsingMatrix::getE() {
	return E;
}

double IsingMatrix::getM() {
	return M;
}

int& IsingMatrix::operator() (int x, int y) {
	if (x > dx) {
		x -= dx;
	} else if (x < 0) {
		x += dx;
	}
	if (y > dy) {
		y -= dy;
	} else if (y < 0) {
		y += dy;
	}
	return matrix[x+dx*y];
}

int* IsingMatrix::operator[] (int i) {
	return &matrix[i];
}
