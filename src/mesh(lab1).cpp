#include "../headers/mesh(lab1).h"

int GlobalData::ReadFromFile()
{
	fstream plik;
	plik.open("data.txt", ios::in);
	if (plik.good() == true)
	{
		while (!plik.eof())
		{
			plik >> local_H;
			plik >> W;
			plik >> nH;
			plik >> nW;
			plik >> npc;
			plik >> k;
			plik >> ro;
			plik >> cp;
			plik >> t0;
		}
		plik.close();
	}
	else {
		cout << "data.txt file is missing, reading default data" << endl;
		system("pause");
		local_H = 0.1;
		W = 0.1;
		nH = 4;
		nW = 4;
		npc = 2;
		k = 25;
		ro = 7800;
		cp = 700;
		t0 = 100;
	}

	return 0;
}

GlobalData::GlobalData()
{
	GlobalData* GB = this;

	GB->ReadFromFile();
	npc = npc * npc;

	dy = GB->local_H / (GB->nH - 1);
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

			if (i == 0 || i == GB->nW-1 || j == 0 || j == GB->nH-1) {
				ND[GB->temp_n].BC = 1;
			}
			else {
			ND[GB->temp_n].BC = 0;
			}
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
		cout << ND[i].x << ", " << ND[i].y << ", " << ND[i].BC << endl;
	}
	cout << endl;
	for (int i = 0; i < GB->nE; i++) {
		cout << Elem[i].ID[0] << ", " << Elem[i].ID[1] << ", " << Elem[i].ID[2] << ", " << Elem[i].ID[3] << endl;
	}
}

int Element::initialize_H_and_C(double xy[2][4], Elem4* e, GlobalData* GB, Node ND[4]) {
	double J[2][2] = { 0.,0.,0.,0. };
	double detJ;

	double wynik_x[4] = { 0.,0.,0.,0. };
	double wynik_y[4] = { 0.,0.,0.,0. };
	double mx[4][4];
	double my[4][4];
	double m_end[4][4];

	for (int x = 0; x < e->npc; x++) {
		//cout << x + 1 << " integration point" << endl;
		//filling J with zeros
		for (int i = 0; i < 2; i++) {
			for (int j = 0; j < 2; j++) {
				J[i][j] = 0;
			}
		}

		fill_J(J, e, xy, x);
		//cout << endl;
		//print_M(J);

		detJ = det_J(J);
		//cout << endl << "det: " << detJ << endl;

		reverse_J(J);
		//cout << "reverse";
		//print_M(J);

		wynik_x[0] = (1. / detJ) * (J[0][0] * e->tab_ksi[x][0] + J[0][1] * e->tab_eta[x][0]);
		wynik_x[1] = (1. / detJ) * (J[0][0] * e->tab_ksi[x][1] + J[0][1] * e->tab_eta[x][1]);
		wynik_x[2] = (1. / detJ) * (J[0][0] * e->tab_ksi[x][2] + J[0][1] * e->tab_eta[x][2]);
		wynik_x[3] = (1. / detJ) * (J[0][0] * e->tab_ksi[x][3] + J[0][1] * e->tab_eta[x][3]);

		/*cout << endl;
		for (int i = 0; i < 4; i++) {
			cout << wynik_x[i] << ", ";
		}
		cout << endl;*/

		for (int i = 0; i < 4; i++) {
			for (int j = 0; j < 4; j++) {
				mx[i][j] = wynik_x[i] * wynik_x[j];
			}
		}

		wynik_y[0] = (1. / detJ) * (J[1][0] * e->tab_ksi[x][0] + J[1][1] * e->tab_eta[x][0]);
		wynik_y[1] = (1. / detJ) * (J[1][0] * e->tab_ksi[x][1] + J[1][1] * e->tab_eta[x][1]);
		wynik_y[2] = (1. / detJ) * (J[1][0] * e->tab_ksi[x][2] + J[1][1] * e->tab_eta[x][2]);
		wynik_y[3] = (1. / detJ) * (J[1][0] * e->tab_ksi[x][3] + J[1][1] * e->tab_eta[x][3]);

		/*cout << endl;
		for (int i = 0; i < 4; i++) {
			cout << wynik_y[i] << ", ";
		}
		cout << endl;*/

		for (int i = 0; i < 4; i++) {
			for (int j = 0; j < 4; j++) {
				my[i][j] = wynik_y[i] * wynik_y[j];
			}
		}

		for (int i = 0; i < 4; i++) {
			for (int j = 0; j < 4; j++) {
				m_end[i][j] = my[i][j] + mx[i][j];
				m_end[i][j] *= detJ * GB->k;
			}
		}

		//print_M(m_end);
		//cout << endl << endl;

		for (int i = 0; i < 4; i++) {
			for (int j = 0; j < 4; j++) {
				local_H[i][j] += m_end[i][j] * e->w1[x] * e->w2[x];
				local_C[i][j] += e->N[x][i] * e->N[x][j] * GB->cp * GB->ro * detJ * e->w1[x] * e->w2[x];
			}
		}
	}

	cout << ND[0].BC << endl;
	cout << ND[1].BC << endl;
	cout << ND[2].BC << endl;
	cout << ND[3].BC << endl;

	double bc_matrix[4][4];
	for (int i = 0; i < 4; i++) {
		for (int j = 0; j < 4; j++) {
			bc_matrix[i][j] = 0;
		}
	}

	if (ND[0].BC == 1 && ND[1].BC == 1) {
		Elem4::Surface* surface = new Elem4::Surface(e->npc, 0);
		for (int h = 0; h < surface->npc; h++) {
			surface->surface_N[0] = surface->N1[h];
			surface->surface_N[1] = surface->N2[h];
			surface->surface_N[2] = 0;
			surface->surface_N[3] = 0;

			for (int i = 0; i < 4; i++) {
				for (int j = 0; j < 4; j++) {
					bc_matrix[i][j] += surface->surface_N[i] * surface->surface_N[j] * 25 * (0.0333333 / 2.);
				}
			}
		}
	}
	

	if (ND[1].BC == 1 && ND[2].BC == 1) {
		Elem4::Surface* surface = new Elem4::Surface(e->npc, 1);
		for (int h = 0; h < surface->npc; h++) {
			surface->surface_N[0] = 0;
			surface->surface_N[1] = surface->N1[h];
			surface->surface_N[2] = surface->N2[h];
			surface->surface_N[3] = 0;

			for (int i = 0; i < 4; i++) {
				for (int j = 0; j < 4; j++) {
					bc_matrix[i][j] += surface->surface_N[i] * surface->surface_N[j] * 25 * (0.0333333 / 2.);
				}
			}
		}
	}

	if (ND[2].BC == 1 && ND[3].BC == 1) {
		Elem4::Surface* surface = new Elem4::Surface(e->npc, 2);
		for (int h = 0; h < surface->npc; h++) {
			surface->surface_N[0] = 0;
			surface->surface_N[1] = 0;
			surface->surface_N[2] = surface->N1[h];
			surface->surface_N[3] = surface->N2[h];

			for (int i = 0; i < 4; i++) {
				for (int j = 0; j < 4; j++) {
					bc_matrix[i][j] += surface->surface_N[i] * surface->surface_N[j] * 25 * (0.0333333 / 2.);
				}
			}
		}
	}

	if (ND[3].BC == 1 && ND[0].BC == 1) {
		Elem4::Surface* surface = new Elem4::Surface(e->npc, 3);
		for (int h = 0; h < surface->npc; h++) {
			surface->surface_N[0] = surface->N2[h];
			surface->surface_N[1] = 0;
			surface->surface_N[2] = 0;
			surface->surface_N[3] = surface->N1[h];

			for (int i = 0; i < 4; i++) {
				for (int j = 0; j < 4; j++) {
					bc_matrix[i][j] += surface->surface_N[i] * surface->surface_N[j] * 25 * (0.0333333 / 2.);
				}
			}
		}
	}
	
	system("pause");

	for (int i = 0; i < 4; i++) {
		for (int j = 0; j < 4; j++) {
			local_H[i][j] += bc_matrix[i][j];
		}
	}

	return 0;
}
Element::Element() {
	//filling H with zeros
	for (int i = 0; i < 4; i++) {
		for (int j = 0; j < 4; j++) {
			local_H[i][j] = 0;
		}
	}

	//filling C with zeros
	for (int i = 0; i < 4; i++) {
		for (int j = 0; j < 4; j++) {
			local_C[i][j] = 0;
		}
	}
}

SoE::SoE(GlobalData* GB){
	//init & filling global_H and global_C with zeros
	global_H = new double* [GB->nN];
	for (int i = 0; i < GB->nN; ++i) {
		global_H[i] = new double[GB->nN];
	}
	global_C = new double* [GB->nN];
	for (int i = 0; i < GB->nN; ++i) {
		global_C[i] = new double[GB->nN];
	}
	for (int k = 0; k < GB->nN; k++) {
		for (int l = 0; l < GB->nN; l++) {
			global_H[k][l] = 0;
			global_C[k][l] = 0;
		}
	}
}
