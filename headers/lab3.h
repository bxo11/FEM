#ifndef LAB3_H
#define LAB3_H

#include <iostream>
#include <cstdlib>

using namespace std;

struct Elem4
{
	int npc;
	double value = 1. / sqrt(3);
	double* ksi;
	double* eta;
	double* w1;
	double* w2;
	double** tab_ksi; // dn/dksi
	double** tab_eta;
	double** N;

	Elem4(int npc);
};

int fill_J(double J[2][2], Elem4* e, double xy[2][4], int i);
double det_J(double J[2][2]);
int reverse_J(double J[2][2]);
void print_M(double J[2][2]);
void print_M(double J[4][4]);

#endif