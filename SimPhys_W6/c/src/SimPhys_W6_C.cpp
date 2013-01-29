//============================================================================
// Name        : SimPhys_W6_C.cpp
// Author      : 
// Version     :
// Copyright   : Your copyright notice
// Description : Hello World in C++, Ansi-style
//============================================================================

#include <iostream>
#include <fstream>
#include <getopt.h>
#include <stdlib.h>
using namespace std;

int main(int argc, char* argv[]) {
	int Ising_L = 4;
	float Ising_J = 1;
	float Ising_H = 0;
	int MC_Seed = 42;
	int MC_Steps = 1000;
	float T_Start = 1;
	float T_Stop = 5;
	float T_StepSize = 0.1;
	int Binning_K = 50;
	char *filename = NULL;

    static struct option long_options[] = {
        {"Ising_L", 1, 0, 0},
        {"Ising_J", 1, 0, 0},
        {"Ising_H", 1, 0, 0},
        {"MC_Seed", 1, 0, 0},
        {"MC_Steps", 1, 0, 0},
        {"T_Start", 1, 0, 0},
        {"T_Stop", 1, 0, 0},
        {"T_StepSize", 1, 0, 0},
        {"Binning_K", 1, 0, 0},
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
				MC_Steps = atof(optarg);
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
			case 8:
				Binning_K = atof(optarg);
				break;
			}
			break;
		case 'o':
			filename = optarg;
			break;
    	}
    }

    std::cout << Ising_L << endl;
    std::cout << Ising_J << endl;
    std::cout << Ising_H << endl;
    std::cout << MC_Seed << endl;
    std::cout << MC_Steps << endl;
    std::cout << T_Start << endl;
    std::cout << T_Stop << endl;
    std::cout << T_StepSize << endl;
    std::cout << Binning_K << endl;

	ofstream myfile(filename);
	if (myfile.is_open())
	{
		myfile << Ising_L<< endl;
		myfile << Ising_J+5<< endl;
		myfile << "3"<< endl;
		myfile.close();
	}
	else cout << "Unable to open file!";

	return 0;
}
