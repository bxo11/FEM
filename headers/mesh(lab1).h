#ifndef LAB1_H
#define LAB1_H

#include <fstream>
#include "lab3.h"
using namespace std;

struct Node
{
	double x, y, t0;
};

struct GlobalData
{
	double dx, dy;
	int temp_n;

	double local_H, W, nH, nW, npc, k, ro, cp, t0; //from file
	double nE, nN;
	int ReadFromFile();
	GlobalData();
};

struct Element
{
	int ID[4];
	double local_H[4][4];
	double local_C[4][4];
	int initialize_H_and_C(double xy[2][4], Elem4* e, GlobalData* GB);

	Element();
};

struct SoE
{
	double** global_H;
	double** global_C;

	SoE(GlobalData* GB);
};

void meshInit(GlobalData* GB, Node* ND, Element* Elem);
void meshPrint(GlobalData* GB, Node* ND, Element* Elem);
#endif