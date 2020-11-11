#include <iostream>
#include <fstream>
#include <cstdlib>
#include <string>
#include "../headers/mesh(lab1).h"
using namespace std;



int main2() {

	/*GlobalData* GB = new GlobalData;

	GB->ReadFromFile();

	dy = GB->H / (GB->nH - 1);
	dx = GB->W / (GB->nW - 1);
	
	
	GB->nE = (GB->nH - 1) * (GB->nW - 1);
	GB->nN = GB->nH * GB->nW;


	Node* ND = new Node[GB->nN];
	temp_n = 0;
	for (int i = 0; i < GB->nW; i++) {
		for (int j = 0; j < GB->nH; j++) {
			ND[temp_n].x = i*dx;
			ND[temp_n].y = j*dy;
			temp_n++;
		}
	}

	Element* Elem = new Element[GB->nE];
	temp_n = 0;
	for (int i = 0; i < GB->nW-1; i++) {
		for (int j = 0; j < GB->nH-1; j++) {
			Elem[temp_n].ID[0] = 1+j+(i*GB->nH);
			Elem[temp_n].ID[1] = Elem[temp_n].ID[0] + GB->nH;
			Elem[temp_n].ID[2] = Elem[temp_n].ID[1] + 1;
			Elem[temp_n].ID[3] = Elem[temp_n].ID[0] + 1;
			temp_n++;
		}
	}

	
	for (int i = 0; i < GB->nN; i++) {
		cout << ND[i].x << ", " <<  ND[i].y << endl;
	}
	cout << endl;
	for (int i = 0; i < GB->nE; i++) {
		cout << Elem[i].ID[0] << ", " << Elem[i].ID[1] << ", " << Elem[i].ID[2] << ", " << Elem[i].ID[3] << endl;
	}*/


	GlobalData* GB = new GlobalData;
	Node* ND = new Node[GB->nN];
	Element* Elem = new Element[GB->nE];
	meshInit(GB, ND, Elem);
	meshPrint(GB, ND, Elem);
	system("pause");
	return 0;
}


