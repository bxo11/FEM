#ifndef LAB1_H
#define LAB1_H

#include <fstream>
#include "part2.h"
using namespace std;

struct Node
{
	double x, y, t0=0;
	int BC;
};

struct GlobalData
{
	double dx, dy;
	int temp_n;

	double H, W, nH, nW, npc, k, ro, cp, t0 , alfa, totoczenia; //from file
	double nE, nN, dTime;
	int ReadFromFile();
	GlobalData();
};

struct Element
{
	int ID[4];
	double local_H[4][4];
	double local_C[4][4];
	double local_P[4];
	int initialize_H_C_P(double xy[2][4], Elem4* e, GlobalData* GB, Node ND[4]);

	Element();
};

struct SoE
{
	double** global_H;
	double** global_C;
	double* global_P;

	SoE(GlobalData* GB);
};

void meshInit(GlobalData* GB, Node* ND, Element* Elem);
void meshPrint(GlobalData* GB, Node* ND, Element* Elem);
#endif