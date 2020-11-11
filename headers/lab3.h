#include <iostream>
#include <fstream>
#include <cstdlib>
#include <string>
using namespace std;

struct Elem4
{
	double value = 1. / sqrt(3);
	double ksi[4] = { -value,value,value,-value };
	double eta[4] = { -value,-value,value,value };
	double tab_ksi[4][4];
	double tab_eta[4][4];

	Elem4() {
		for (int i = 0; i < 4; i++) {
			tab_ksi[i][0] = -1. / 4 * (1 - eta[i]);
			tab_ksi[i][1] = 1. / 4 * (1 - eta[i]);
			tab_ksi[i][2] = 1. / 4 * (1 + eta[i]);
			tab_ksi[i][3] = -1. / 4 * (1 + eta[i]);

			tab_eta[i][0] = -1. / 4 * (1 - ksi[i]);
			tab_eta[i][1] = -1. / 4 * (1 + ksi[i]);
			tab_eta[i][2] = 1. / 4 * (1 + ksi[i]);
			tab_eta[i][3] = 1. / 4 * (1 - ksi[i]);
		}
	}
};

int fill_J(double J[2][2], Elem4* e, double xy[2][4], int i);
double det_J(double J[2][2]);
int reverse_J(double J[2][2]);
void print_M(double J[2][2]);
void print_M(double J[4][4]);

int fill_J(double J[2][2], Elem4* e, double xy[2][4], int i) {
	for (int k = 0; k < 4; k++) {
		J[0][0] += e->tab_ksi[i][k] * xy[0][k];
		J[0][1] += e->tab_ksi[i][k] * xy[1][k];
		J[1][0] += e->tab_eta[i][k] * xy[0][k];
		J[1][1] += e->tab_eta[i][k] * xy[1][k];
	}
	return 0;
}

inline double det_J(double J[2][2])
{
	return (J[0][0] * J[1][1]) - (J[0][1] * J[1][0]);
}

inline int reverse_J(double J[2][2])
{
	swap(J[0][0], J[1][1]);
	J[0][1] *= -1.;
	J[1][0] *= -1.;
	return 0;
}

inline void print_M(double J[2][2])
{
	cout << "Jacobi";
	cout << endl << "[ " << J[0][0] << ", " << J[0][1];
	cout << endl << "  " << J[1][0] << ", " << J[1][1] << " ]";
}

inline void print_M(double M[4][4])
{
	for (int i = 0; i < 4; i++) {
		for (int j = 0; j < 4; j++) {
			if (j == 0) cout << endl;
			cout << M[i][j] << ", ";
		}
	}
}
