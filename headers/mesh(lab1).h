#ifndef LAB1_H
#define LAB1_H

#include <fstream>
#include "lab3.h"
using namespace std;

struct Node
{
	double x, y;
};

struct GlobalData
{
	double dx, dy;
	int temp_n;

	double H, W, nH, nW, npc, k; //from file
	double nE, nN;
	int ReadFromFile();
	GlobalData();
};

struct Element
{
	int ID[4];
	double H[4][4];
	int initialize_H(double xy[2][4], Elem4* e, GlobalData* GB);

	Element();
};

void meshInit(GlobalData* GB, Node* ND, Element* Elem);
void meshPrint(GlobalData* GB, Node* ND, Element* Elem);
#endif