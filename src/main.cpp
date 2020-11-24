#include "../headers/lab3.h"
#include "../headers/mesh(lab1).h"

using namespace std;

int main()
{
	GlobalData* GB = new GlobalData;
	Node* ND = new Node[GB->nN];
	Element* Elem = new Element[GB->nE];
	meshInit(GB, ND, Elem);

	Elem4* e = new Elem4(GB->npc);
	SoE* soe = new SoE(GB);

	for (int i = 0; i < GB->nE; i++) {
		double xy[2][4] = { ND[Elem[i].ID[0]].x,ND[Elem[i].ID[1]].x,ND[Elem[i].ID[2]].x,ND[Elem[i].ID[3]].x,
							ND[Elem[i].ID[0]].y,ND[Elem[i].ID[1]].y,ND[Elem[i].ID[2]].y,ND[Elem[i].ID[3]].y, };

		Elem[i].initialize_H_and_C(xy, e,GB);

		for (int k = 0; k < 4; k++) {
			for (int l = 0; l < 4; l++) {
				soe->global_H[Elem[i].ID[k]][Elem[i].ID[l]] += Elem[i].local_H[k][l];
				soe->global_C[Elem[i].ID[k]][Elem[i].ID[l]] += Elem[i].local_C[k][l];
			}
		}
	}

	//system("cls");
	meshPrint(GB, ND, Elem);

	cout << endl;
	for (int k = 0; k < GB->nN; k++) {
		for (int l = 0; l < GB->nN; l++) {
			cout << soe->global_H[k][l] << ", ";
		}
		cout << endl;
	}

	cout << endl;
	for (int k = 0; k < GB->nN; k++) {
		for (int l = 0; l < GB->nN; l++) {
			cout << soe->global_C[k][l] << ", ";
		}
		cout << endl;
	}

	return 0;
}