#include "../headers/lab3.h"
#include "../headers/mesh(lab1).h"

Elem4::Elem4(int a_npc) {
	this->npc = a_npc;

	tab_ksi = new double* [npc];
	for (int i = 0; i < npc; ++i) {
		tab_ksi[i] = new double[4];
	}

	tab_eta = new double* [npc];
	for (int i = 0; i < npc; ++i) {
		tab_eta[i] = new double[4];
	}

	N = new double* [npc];
	for (int i = 0; i < npc; ++i) {
		N[i] = new double[4];
	}

	for (int k = 0; k < npc; k++) {
		for (int l = 0; l < 4; l++) {
			tab_ksi[k][l] = 0;
			tab_eta[k][l] = 0;
			N[k][l] = 0;
		}
	}

	double value;
	switch (npc)
	{
	case 4:

		value = 1. / sqrt(3);
		ksi = new double[4]{ -value,value,value,-value };
		eta = new double[4]{ -value,-value,value,value };
		w1 = new double[4]{ 1,1,1,1 };
		w2 = new double[4]{ 1,1,1,1 };
		break;

	case 9:

		value = sqrt(3. / 5.);
		ksi = new double[9]{ -value,0.,value,value,0.,-value,value,0.,-value };
		eta = new double[9]{ -value,-value,-value,0.,0.,0.,value,value,value };
		w1 = new double[9]{ 5. / 9, 5. / 9, 5. / 9, 8. / 9, 8. / 9, 8. / 9, 5. / 9, 5. / 9, 5. / 9 };
		w2 = new double[9]{ 5. / 9, 8./9, 5. / 9, 5. / 9, 8. / 9, 5. / 9, 5. / 9, 8. / 9, 5. / 9 };
		break;

	default:
		cout << "Error - integration points" << endl;
		system("pause");
		break;
	}

	for (int i = 0; i < npc; i++) {
		tab_ksi[i][0] = -1. / 4 * (1 - eta[i]);
		tab_ksi[i][1] = 1. / 4 * (1 - eta[i]);
		tab_ksi[i][2] = 1. / 4 * (1 + eta[i]);
		tab_ksi[i][3] = -1. / 4 * (1 + eta[i]);

		tab_eta[i][0] = -1. / 4 * (1 - ksi[i]);
		tab_eta[i][1] = -1. / 4 * (1 + ksi[i]);
		tab_eta[i][2] = 1. / 4 * (1 + ksi[i]);
		tab_eta[i][3] = 1. / 4 * (1 - ksi[i]);

		N[i][0] = 1. / 4 * (1 - ksi[i]) * (1 - eta[i]);
		N[i][1] = 1. / 4 * (1 + ksi[i]) * (1 - eta[i]);
		N[i][2] = 1. / 4 * (1 + ksi[i]) * (1 + eta[i]);
		N[i][3] = 1. / 4 * (1 - ksi[i]) * (1 + eta[i]);
	}

	//for (int k = 0; k < npc; k++) {
	//	for (int l = 0; l < 4; l++) {
	//		cout << tab_ksi[k][l] << ", ";
	//		//tab_eta[k][l] = 0;
	//	}
	//	cout << endl;
	//}
}

Elem4::Surface::Surface(int a_npc, int x){
	this->npc = sqrt(a_npc);

	switch (x)
	{
	case 0:
		switch (npc)
		{
		case 2:
			pcKsi = new double[npc] {-1. / sqrt(3), 1. / sqrt(3)};
			pcEta = new double[npc] {-1., -1.};
			wpc = new double[npc] {1., 1.};
			break;
		case 3:
			pcKsi = new double[npc] {-sqrt(3. / 5), 0., sqrt(3. / 5)};
			pcEta = new double[npc] {-1., -1., -1.};
			wpc = new double[npc] {5. / 9, 8. / 9, 5. / 9};
			break;
		default:
			break;
		}

		N1 = new double[npc];
		N2 = new double[npc];
		for (int i = 0; i < npc; i++) {
			N1[i] = 1. / 4 * (1 - pcKsi[i]) * (1 - pcEta[i]);
			N2[i] = 1. / 4 * (1 + pcKsi[i]) * (1 - pcEta[i]);
		}

		break;
	case 1:
		switch (npc)
		{
		case 2:
			pcKsi = new double[npc] {1., 1.};
			pcEta = new double[npc] {-1. / sqrt(3), 1. / sqrt(3)};
			wpc = new double[npc] {1., 1.};
			break;
		case 3:
			pcKsi = new double[npc] {1., 1., 1.};
			pcEta = new double[npc] {-sqrt(3. / 5), 0., sqrt(3. / 5)};
			wpc = new double[npc] {5. / 9, 8. / 9, 5. / 9};
			break;
		default:
			break;
		}

		N1 = new double[npc];
		N2 = new double[npc];
		for (int i = 0; i < npc; i++) {
			N1[i] = 1. / 4 * (1 + pcKsi[i]) * (1 - pcEta[i]);
			N2[i] = 1. / 4 * (1 + pcKsi[i]) * (1 + pcEta[i]);
		}

		break;
	case 2:
		switch (npc)
		{
		case 2:
			pcKsi = new double[npc] {1. / sqrt(3), -1. / sqrt(3)};
			pcEta = new double[npc] {1., 1.};
			wpc = new double[npc] {1., 1.};
			break;
		case 3:
			pcKsi = new double[npc] {sqrt(3. / 5), 0., -sqrt(3. / 5)};
			pcEta = new double[npc] {1., 1., 1.};
			wpc = new double[npc] {5. / 9, 8. / 9, 5. / 9};
			break;
		default:
			break;
		}

		N1 = new double[npc];
		N2 = new double[npc];
		for (int i = 0; i < npc; i++) {
			N1[i] = 1. / 4 * (1 + pcKsi[i]) * (1 + pcEta[i]);
			N2[i] = 1. / 4 * (1 - pcKsi[i]) * (1 + pcEta[i]);
		}

		break;
	case 3:
		switch (npc)
		{
		case 2:
			pcKsi = new double[npc] {-1., -1.};
			pcEta = new double[npc] {1. / sqrt(3), -1. / sqrt(3)};
			wpc = new double[npc] {1., 1.};
			break;
		case 3:
			pcKsi = new double[npc] {-1., -1., -1.};
			pcEta = new double[npc] {sqrt(3. / 5), 0., -sqrt(3. / 5)};
			wpc = new double[npc] {5. / 9, 8. / 9, 5. / 9};
			break;
		default:
			break;
		}

		N1 = new double[npc];
		N2 = new double[npc];
		for (int i = 0; i < npc; i++) {
			N1[i] = 1. / 4 * (1 - pcKsi[i]) * (1 + pcEta[i]);
			N2[i] = 1. / 4 * (1 - pcKsi[i]) * (1 - pcEta[i]);
		}

		break;
	default:
		break;
	}
	/*
	switch (this->npc)
	{
	case 2:
		switch (x)
		{
		case 0:
			pcKsi = new double[npc] {-1. / sqrt(3), 1. / sqrt(3)};
			pcEta = new double[npc] {-1., -1.};
			wpc = new double[npc] {1., 1.};
			N1 = new double[npc];
			N2 = new double[npc];
			for (int i = 0; i < npc; i++) {
				N1[i] = 1. / 4 * (1 - pcKsi[i]) * (1 - pcEta[i]);
				N2[i] = 1. / 4 * (1 + pcKsi[i]) * (1 - pcEta[i]);
			}
			//N1 = new double[npc] {1. / 4 * (1 - pcKsi[0]) * (1 - pcEta[0]), 1. / 4 * (1 - pcKsi[1]) * (1 - pcEta[1]) };
			//N2 = new double[npc] {1. / 4 * (1 + pcKsi[0]) * (1 - pcEta[0]), 1. / 4 * (1 + pcKsi[1]) * (1 - pcEta[1]) };
			break;

		case 1:
			pcKsi = new double[npc] {1., 1.};
			pcEta = new double[npc] {-1. / sqrt(3), 1. / sqrt(3)};
			wpc = new double[npc] {1., 1.};
			N1 = new double[npc];
			N2 = new double[npc];
			for (int i = 0; i < npc; i++) {
				N1[i] = 1. / 4 * (1 + pcKsi[i]) * (1 - pcEta[i]);
				N2[i] = 1. / 4 * (1 + pcKsi[i]) * (1 + pcEta[i]);
			}
			//N1 = new double[npc] {1. / 4 * (1 + pcKsi[0]) * (1 - pcEta[0]), 1. / 4 * (1 + pcKsi[1]) * (1 - pcEta[1]) };
			//N2 = new double[npc] {1. / 4 * (1 + pcKsi[0]) * (1 + pcEta[0]), 1. / 4 * (1 + pcKsi[1]) * (1 + pcEta[1]) };
			break;

		case 2:
			pcKsi = new double[npc] {1. / sqrt(3), -1. / sqrt(3)};
			pcEta = new double[npc] {1., 1.};
			wpc = new double[npc] {1., 1.};
			N1 = new double[npc];
			N2 = new double[npc];
			for (int i = 0; i < npc; i++) {
				N1[i] = 1. / 4 * (1 + pcKsi[i]) * (1 + pcEta[i]);
				N2[i] = 1. / 4 * (1 - pcKsi[i]) * (1 + pcEta[i]);
			}
			//N1 = new double[npc] {1. / 4 * (1 + pcKsi[0]) * (1 + pcEta[0]), 1. / 4 * (1 + pcKsi[1]) * (1 + pcEta[1]) };
			//N2 = new double[npc] {1. / 4 * (1 - pcKsi[0]) * (1 + pcEta[0]), 1. / 4 * (1 - pcKsi[1]) * (1 + pcEta[1]) };
			break;

		case 3:
			pcKsi = new double[npc] {-1., -1.};
			pcEta = new double[npc] {1. / sqrt(3), -1. / sqrt(3)};
			wpc = new double[npc] {1., 1.};
			N1 = new double[npc];
			N2 = new double[npc];
			for (int i = 0; i < npc; i++) {
				N1[i] = 1. / 4 * (1 - pcKsi[i]) * (1 + pcEta[i]);
				N2[i] = 1. / 4 * (1 - pcKsi[i]) * (1 - pcEta[i]);
			}
			//N1 = new double[npc] {1. / 4 * (1 - pcKsi[0]) * (1 + pcEta[0]), 1. / 4 * (1 - pcKsi[1]) * (1 + pcEta[1]) };
			//N2 = new double[npc] {1. / 4 * (1 - pcKsi[0]) * (1 - pcEta[0]), 1. / 4 * (1 - pcKsi[1]) * (1 - pcEta[1]) };
			break;

		default:
			break;
		}
		break;

	case 3:

		break;

	default:
		break;
	}
	*/
}

int fill_J(double J[2][2], Elem4* e, double xy[2][4], int i) {
	for (int k = 0; k < 4; k++) {
		J[0][0] += e->tab_ksi[i][k] * xy[0][k];
		J[0][1] += e->tab_ksi[i][k] * xy[1][k];
		J[1][0] += e->tab_eta[i][k] * xy[0][k];
		J[1][1] += e->tab_eta[i][k] * xy[1][k];
	}
	return 0;
}

double det_J(double J[2][2])
{
	return (J[0][0] * J[1][1]) - (J[0][1] * J[1][0]);
}

int reverse_J(double J[2][2])
{
	swap(J[0][0], J[1][1]);
	J[0][1] *= -1.;
	J[1][0] *= -1.;
	return 0;
}

void print_M(double J[2][2])
{
	cout << "Jacobi";
	cout << endl << "[ " << J[0][0] << ", " << J[0][1];
	cout << endl << "  " << J[1][0] << ", " << J[1][1] << " ]";
}

void print_M(double M[4][4])
{
	for (int i = 0; i < 4; i++) {
		for (int j = 0; j < 4; j++) {
			if (j == 0) cout << endl;
			cout << M[i][j] << ", ";
		}
	}
}


