﻿#include "../headers/part2.h"
#include "../headers/part1.h"

#define iterations 20
#define tau 1

using namespace std;

int main()
{
	GlobalData* GB = new GlobalData;
	Node* ND = new Node[GB->nN];
	Element* Elem = new Element[GB->nE];
	meshInit(GB, ND, Elem);

	Elem4* e = new Elem4(GB->npc);
	SoE* soe = new SoE(GB);
	double* tx = new double[GB->nN];

	double Asr = GB->k / (GB->cp * GB->ro);
	GB->dTime = pow(GB->W / GB->nW, 2.) / (0.5 * Asr);

	for (int z = 0; z < iterations; z++) {

		for (int k = 0; k < GB->nN; k++) {
			for (int l = 0; l < GB->nN; l++) {
				soe->global_H[k][l] = 0;
				soe->global_C[k][l] = 0;
			}
		}

		//filling global_P with zeros
		for (int i = 0; i < GB->nN; ++i) {
			soe->global_P[i] = 0;
		}

		for (int i = 0; i < GB->nE; i++) {
			double xy[2][4] = { ND[Elem[i].ID[0]].x,ND[Elem[i].ID[1]].x,ND[Elem[i].ID[2]].x,ND[Elem[i].ID[3]].x,
								ND[Elem[i].ID[0]].y,ND[Elem[i].ID[1]].y,ND[Elem[i].ID[2]].y,ND[Elem[i].ID[3]].y, };

			Node ND_tab[4] = { ND[Elem[i].ID[0]] ,ND[Elem[i].ID[1]] ,ND[Elem[i].ID[2]] ,ND[Elem[i].ID[3]] };
			Elem[i].initialize_H_C_P(xy, e, GB, ND_tab);
			//print_M(Elem[i].local_C);
			//cout << endl;
			for (int k = 0; k < 4; k++) {
				for (int l = 0; l < 4; l++) {
					soe->global_H[Elem[i].ID[k]][Elem[i].ID[l]] += Elem[i].local_H[k][l];
					soe->global_C[Elem[i].ID[k]][Elem[i].ID[l]] += Elem[i].local_C[k][l];
				}
			}

			for (int k = 0; k < 4; k++) {
				soe->global_P[Elem[i].ID[k]] += Elem[i].local_P[k];
			}
		}

		//system("cls");
		//meshPrint(GB, ND, Elem);

		/*cout << endl;
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
		}*/

		for (int k = 0; k < GB->nN; k++) {
			for (int l = 0; l < GB->nN; l++) {
				soe->global_H[k][l] += soe->global_C[k][l] / tau;
			}
		}

		//cout << endl;
		//for (int k = 0; k < GB->nN; k++) {
		//	for (int l = 0; l < GB->nN; l++) {
		//		cout << soe->global_H[k][l] << ", ";
		//	}
		//	cout << endl;
		//}

		double* temp_tab = new double[GB->nN];
		for (int t = 0; t < GB->nN; t++) {
			temp_tab[t] = 0;
		}

		for (int k = 0; k < GB->nN; k++) {
			soe->global_P[k] *= -1;
		}
		for (int k = 0; k < GB->nN; k++) {

			for (int l = 0; l < GB->nN; l++) {
				temp_tab[k] += ND[l].t0 * soe->global_C[k][l] / tau;
			}
			soe->global_P[k] += temp_tab[k];
		}

		/*cout << endl;
		for (int k = 0; k < GB->nN; k++) {
			cout << soe->global_P[k] << ", ";
		}*/

		tx = Gauss_elimination(GB->nN, soe->global_H, soe->global_P);

		double min=1e9;
		double max=-1e9;
		for (int i = 0; i < GB->nN; i++) {
			if (tx[i] > max) {
				max = tx[i];
			}
			if (tx[i] < min) {
				min=tx[i];
			}
		}
		cout <<"Iteration "<<z<< " - Min & Max: " << min << ", " << max << endl;

		for (int i = 0; i < GB->nN; i++) {
			ND[i].t0 = tx[i];
		}
	}
	return 0;
}