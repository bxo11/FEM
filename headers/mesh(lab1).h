#include <iostream>
#include <fstream>
#include <cstdlib>
#include <string>
#include "lab3.h"
using namespace std;

struct Node
{
	double x, y;
};

struct Element
{
	int ID[4];
	Node* ND;
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



