#ifndef LAB3_H
#define LAB3_H

#include <iostream>
#include <cstdlib>

using namespace std;

struct Elem4
{
	double value = 1. / sqrt(3);
	double ksi[4] = { -value,value,value,-value };
	double eta[4] = { -value,-value,value,value };
	double tab_ksi[4][4];
	double tab_eta[4][4];

	Elem4();
};

int fill_J(double J[2][2], Elem4* e, double xy[2][4], int i);
double det_J(double J[2][2]);
int reverse_J(double J[2][2]);
void print_M(double J[2][2]);
void print_M(double J[4][4]);

#endif