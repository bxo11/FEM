#include <fstream>
#include "lab3.h"
using namespace std;

struct Node
{
	double x, y;
};

struct Element
{
	int ID[4];
	double H[4][4];
	int initialize_H(double xy[2][4], Elem4* e);

	Element() {
		//zerowanie H
		for (int i = 0; i < 4; i++) {
			for (int j = 0; j < 4; j++) {
				H[i][j] = 0;
			}
		}
	}
};

struct GlobalData
{
	double dx, dy;
	int temp_n;

	double H, W, nH, nW; //from file
	double nE, nN;
	int ReadFromFile();
	GlobalData();
};

void meshInit(GlobalData* GB, Node* ND, Element* Elem);
void meshPrint(GlobalData* GB, Node* ND, Element* Elem);

int GlobalData::ReadFromFile()
{
	fstream plik;
	plik.open("data.txt", ios::in);
	if (plik.good() == true)
	{
		while (!plik.eof())
		{
			plik >> H;
			plik >> W;
			plik >> nH;
			plik >> nW;
		}
		plik.close();
	}

	return 0;
}

GlobalData::GlobalData()
{
	GlobalData* GB = this;

	GB->ReadFromFile();

	dy = GB->H / (GB->nH - 1);
	dx = GB->W / (GB->nW - 1);

	GB->nE = (GB->nH - 1) * (GB->nW - 1);
	GB->nN = GB->nH * GB->nW;
}

void meshInit(GlobalData* GB, Node* ND, Element* Elem)
{
	GB->temp_n = 0;
	for (int i = 0; i < GB->nW; i++) {
		for (int j = 0; j < GB->nH; j++) {
			ND[GB->temp_n].x = i * GB->dx;
			ND[GB->temp_n].y = j * GB->dy;
			GB->temp_n++;
		}
	}

	GB->temp_n = 0;
	for (int i = 0; i < GB->nW - 1; i++) {
		for (int j = 0; j < GB->nH - 1; j++) {
			Elem[GB->temp_n].ID[0] = 0 + j + (i * GB->nH);
			Elem[GB->temp_n].ID[1] = Elem[GB->temp_n].ID[0] + GB->nH;
			Elem[GB->temp_n].ID[2] = Elem[GB->temp_n].ID[1] + 1;
			Elem[GB->temp_n].ID[3] = Elem[GB->temp_n].ID[0] + 1;
			GB->temp_n++;
		}
	}
}

void meshPrint(GlobalData* GB, Node* ND, Element* Elem)
{
	for (int i = 0; i < GB->nN; i++) {
		cout << ND[i].x << ", " << ND[i].y << endl;
	}
	cout << endl;
	for (int i = 0; i < GB->nE; i++) {
		cout << Elem[i].ID[0] << ", " << Elem[i].ID[1] << ", " << Elem[i].ID[2] << ", " << Elem[i].ID[3] << endl;
	}
}

int Element::initialize_H(double xy[2][4], Elem4* e) {
	double J[2][2] = { 0.,0.,0.,0. };
	double detJ;

	double wynik_x[4] = { 0.,0.,0.,0. };
	double wynik_y[4] = { 0.,0.,0.,0. };
	double mx[4][4];
	double my[4][4];
	double m_end[4][4];

	for (int x = 0; x < 4; x++) {
		cout << x + 1 << " punkt calkowania" << endl;
		//zerowanie J1
		for (int i = 0; i < 2; i++) {
			for (int j = 0; j < 2; j++) {
				J[i][j] = 0;
			}
		}

		fill_J(J, e, xy, x);
		print_M(J);

		detJ = det_J(J);
		cout << endl << "det: " << detJ << endl;

		reverse_J(J);
		cout << "reverse";
		print_M(J);

		wynik_x[0] = (1. / detJ) * (J[0][0] * e->tab_ksi[x][0] + J[0][1] * e->tab_eta[x][0]);
		wynik_x[1] = (1. / detJ) * (J[0][0] * e->tab_ksi[x][1] + J[0][1] * e->tab_eta[x][1]);
		wynik_x[2] = (1. / detJ) * (J[0][0] * e->tab_ksi[x][2] + J[0][1] * e->tab_eta[x][2]);
		wynik_x[3] = (1. / detJ) * (J[0][0] * e->tab_ksi[x][3] + J[0][1] * e->tab_eta[x][3]);

		cout << endl;
		for (int i = 0; i < 4; i++) {
			cout << wynik_x[i] << ", ";
		}
		cout << endl;

		for (int i = 0; i < 4; i++) {
			for (int j = 0; j < 4; j++) {
				mx[i][j] = wynik_x[i] * wynik_x[j];
			}
		}

		wynik_y[0] = (1. / detJ) * (J[1][0] * e->tab_ksi[x][0] + J[1][1] * e->tab_eta[x][0]);
		wynik_y[1] = (1. / detJ) * (J[1][0] * e->tab_ksi[x][1] + J[1][1] * e->tab_eta[x][1]);
		wynik_y[2] = (1. / detJ) * (J[1][0] * e->tab_ksi[x][2] + J[1][1] * e->tab_eta[x][2]);
		wynik_y[3] = (1. / detJ) * (J[1][0] * e->tab_ksi[x][3] + J[1][1] * e->tab_eta[x][3]);

		cout << endl;
		for (int i = 0; i < 4; i++) {
			cout << wynik_y[i] << ", ";
		}
		cout << endl;

		for (int i = 0; i < 4; i++) {
			for (int j = 0; j < 4; j++) {
				my[i][j] = wynik_y[i] * wynik_y[j];
			}
		}

		for (int i = 0; i < 4; i++) {
			for (int j = 0; j < 4; j++) {
				m_end[i][j] = my[i][j] + mx[i][j];
				m_end[i][j] *= 25 * detJ; //k(wspolczynnik przewodzenia ciepla)=30
			}
		}

		print_M(m_end);
		cout << endl << endl;

		for (int i = 0; i < 4; i++) {
			for (int j = 0; j < 4; j++) {
				H[i][j] += m_end[i][j];
			}
		}
	}

	return 0;
}
