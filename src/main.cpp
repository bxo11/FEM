#include "../headers/lab3.h"
#include "../headers/mesh(lab1).h"

using namespace std;
double** global_H;

int main()
{
	GlobalData* GB = new GlobalData;
	Node* ND = new Node[GB->nN];
	Element* Elem = new Element[GB->nE];
	meshInit(GB, ND, Elem);
	//meshPrint(GB, ND, Elem);

	Elem4* e = new Elem4;

	global_H = new double* [GB->nN];
	for (int i = 0; i < GB->nN; ++i) {
		global_H[i] = new double[GB->nN];
	}

	for (int k = 0; k < GB->nN; k++) {
		for (int l = 0; l < GB->nN; l++) {
			global_H[k][l] = 0;
		}
	}

	for (int i = 0; i < GB->nE; i++) {
		//double xy[2][4] = { 0.,4.,4.,0.,
		//0., 0., 6., 6.

		double xy[2][4] = { ND[Elem[i].ID[0]].x,ND[Elem[i].ID[1]].x,ND[Elem[i].ID[2]].x,ND[Elem[i].ID[3]].x,
							ND[Elem[i].ID[0]].y,ND[Elem[i].ID[1]].y,ND[Elem[i].ID[2]].y,ND[Elem[i].ID[3]].y, };

		Elem[i].initialize_H(xy, e);

		for (int k = 0; k < 4; k++) {
			for (int l = 0; l < 4; l++) {
				global_H[Elem[i].ID[k]][Elem[i].ID[l]] += Elem[i].H[k][l];
			}
		}
	}

	//system("cls");
	meshPrint(GB, ND, Elem);

	for (int k = 0; k < GB->nN; k++) {
		for (int l = 0; l < GB->nN; l++) {
			cout << global_H[k][l] << ", ";
		}
		cout << endl;
	}

	return 0;
}