#include <iostream>
#include <fstream>
#include <cmath>
#include <getopt.h>
#include <stdlib.h>
#include <gsl/gsl_rng.h> //libgsl0-dev
#include <omp.h>
#include <vector>
#include "IsingMatrix.h"
using namespace std;

struct Data {
	Data(int s) : E(s), M(s), T(0), acceptance(0) { }
	vector<double> E, M;
	double T, acceptance;
};

int main(int argc, char* argv[]) {
	//====================
	//=== SETUP ==========
	//====================

	// Parameter
	int Ising_L = 4;
	float Ising_J = 1;
	float Ising_H = 0;
	int MC_Seed = 42;
	int  MC_Sweeps = 100;
	float T_Start = 1;
	float T_Stop = 5;
	float T_StepSize = 0.1;
	char *filename = NULL;

    static struct option long_options[] = {
        {"Ising_L", 1, 0, 0},
        {"Ising_J", 1, 0, 0},
        {"Ising_H", 1, 0, 0},
        {"MC_Seed", 1, 0, 0},
        {"MC_Sweeps", 1, 0, 0},
        {"T_Start", 1, 0, 0},
        {"T_Stop", 1, 0, 0},
        {"T_StepSize", 1, 0, 0},
        {"out", 1, 0, 'o'},
        {NULL, 0, NULL, 0}
    };

    int c, option_index = 0;
    while ((c = getopt_long(argc, argv, "o:", long_options, &option_index)) != -1) {
    	switch (c) {
		case 0:
			switch (option_index) {
			case 0:
				Ising_L = atof(optarg);
				break;
			case 1:
				Ising_J = atof(optarg);
				break;
			case 2:
				Ising_H = atof(optarg);
				break;
			case 3:
				MC_Seed = atof(optarg);
				break;
			case 4:
				MC_Sweeps = atof(optarg);
				break;
			case 5:
				T_Start = atof(optarg);
				break;
			case 6:
				T_Stop = atof(optarg);
				break;
			case 7:
				T_StepSize = atof(optarg);
				break;
			}
			break;
		case 'o':
			filename = optarg;
			break;
    	}
    }

    int V = Ising_L*Ising_L;
    int T_Steps = (T_Stop-T_Start)/T_StepSize+1;

    gsl_rng *rnd;
    rnd = gsl_rng_alloc(gsl_rng_taus);
    gsl_rng_set(rnd, MC_Seed);

    vector<Data> data(T_Steps, Data(MC_Sweeps));
    Data *ptr = &data[0];

    int threads = 4;
    omp_set_num_threads(threads);

    //====================
	//=== CALCULATION ====
	//====================

	#pragma omp parallel for
	for (int l = 0; l < threads; l++) {
		IsingMatrix matrix(Ising_L,Ising_L,MC_Seed,Ising_J,Ising_H);
		matrix.init();
		for (int m = l; m < T_Steps; m+=threads) {
			int n = T_Steps-1-m;
			double T = T_StepSize*n + T_Start;
			double actrate = 0;

			for (int i = 0; i < MC_Sweeps; i++) {
				double sumE = 0;
				double sumM = 0;

				for (int j = 0; j < V; j++) {
					double Delta = matrix.getE();
					matrix.flip();
					Delta -= matrix.getE();

					double r = gsl_rng_uniform(rnd);
					if (r < min(1.0,exp(Delta/T))) {
						actrate++;
					} else {
						matrix.unflip();
					}

					sumE += matrix.getE();
					sumM += matrix.getM();
				}
				ptr[n].E[i] = sumE/V;
				ptr[n].M[i] = abs(sumM/V);
			}

			ptr[n].T = T;
			ptr[n].acceptance = actrate/MC_Sweeps;
		}
	}

	gsl_rng_free(rnd);

	//====================
	//=== WRITING ========
	//====================

	ofstream myfile(filename);
	if (myfile.is_open())
	{
		for (int n = 0; n < T_Steps; n++) {
			myfile << ptr[n].T << " "<< ptr[n].acceptance << " ";
			for (int i = 0; i < MC_Sweeps; i++) {
				myfile << ptr[n].E[i] << " ";
			}
			myfile << endl;
			myfile << ptr[n].T << " "<< ptr[n].acceptance << " ";
			for (int i = 0; i < MC_Sweeps; i++) {
				myfile << ptr[n].M[i] << " ";
			}
			myfile << endl;
		}
		myfile.close();
	}
	else cout << "Unable to open file!";

	return 0;
}
