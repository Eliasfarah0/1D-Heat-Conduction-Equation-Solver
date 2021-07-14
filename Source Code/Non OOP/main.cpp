/***********************************************************************\
#                                                                       #
#                C R A N F I E L D   U N I V E R S I T Y                #
#                          2 0 1 9  /  2 0 2 0                          #
#                                                                       #
#               MSc in Aerospace Computational Engineering              #
#                                                                       #
#                   1D HEAT CONDUCTION EQUATION SOLVER                  #
#                                                                       #
#-----------------------------------------------------------------------#
#                                                                       #
#   Main Contributors:                                                  #
#       Sevan Retif            (Email: Sevan.Retif@cranfield.ac.uk)     #
#       Elias Farah            (Email: E.Farah@cranfield.ac.uk)         #
#                                                                       #
\***********************************************************************/


#include <iostream>
#include <iomanip>
#include <cmath>
#include <fstream>
#include "vector.h"

using namespace std;

const double pi = 3.1415;

Vector Exact_Solution                  (double exact[], double delta_x, double delta_t, int nodes, double output_time, double length, double Diff, double T_sur, double T_init);
Vector DuFort_Frankel_Explicit_Scheme  (double** numerical, double delta_x, double delta_t, double output_time, int nodes, double r, double T_sur, double T_init);
Vector Richardson_Explicit_Scheme      (double** numerical, double delta_x, double delta_t, double output_time, int nodes, double r, double T_sur, double T_init);
Vector Laasonen_Simple_Implicit_Scheme (double** numerical, double delta_x, double delta_t, double output_time, int nodes, double r, double T_sur, double T_init);
Vector Crank_Nicholson_Implicit_Scheme (double** numerical, double delta_x, double delta_t, double output_time, int nodes, double r, double T_sur, double T_init);
int    TDMA_Solver                     (double* lower_diag, double* main_diag, double* upper_diag, double* b, int nodes);
void   Print_Solution                  (Vector v1, Vector v2, Vector v3, Vector v4, Vector v5, int nodes, double delta_x);


int main() {
	double length = 31.0;
	double delta_x	= 0.05;
	double Diff = 93.0;
	double T_sur = 149.0;
	double T_init = 38.0;
	double delta_t = 0.01;
	double output_time = 0.5;

	int nodes = length / delta_x + 1;
	double r = (Diff * delta_t) / (pow(delta_x, 2));

	double* real;
	real = new double[nodes];
	
	double** numerical = new double*[nodes];
	for (int i = 0; i < nodes; i++) numerical[i] = new double[int((output_time / delta_t) + 1)];

	Vector v1, v2, v3, v4, v5;

	v1 = Exact_Solution(real, delta_x, delta_t, nodes, output_time, length, Diff, T_sur, T_init);
	v2 = DuFort_Frankel_Explicit_Scheme(numerical, delta_x, delta_t, output_time, nodes, r, T_sur, T_init);
	v3 = Richardson_Explicit_Scheme(numerical, delta_x, delta_t, output_time, nodes, r, T_sur, T_init);
	v4 = Laasonen_Simple_Implicit_Scheme(numerical, delta_x, delta_t, output_time, nodes, r, T_sur, T_init);
	v5 = Crank_Nicholson_Implicit_Scheme(numerical, delta_x, delta_t, output_time, nodes, r, T_sur, T_init);
	Print_Solution(v1, v2, v3, v4, v5, nodes, delta_x);

	for (int i = 0; i < nodes; i++) delete[] numerical[i];
	delete[] numerical;

	cout << "Computation Completed!" << endl;
	system("PAUSE");;
	return 0;
}


Vector Exact_Solution(double exact[], double delta_x, double delta_t, int nodes, double output_time, double length, double Diff, double T_sur, double T_init) {
	Vector v1;
	double series;
	int acc = 100;

	for (int i = 0; i < nodes; i++) {
		series = 0;
		for (int m = 1; m <= acc; m++) {
			series += (exp(-Diff * output_time * pow((m * pi) / length, 2)) * ((1.0 - pow(-1.0, m)) / (m * pi)) * sin((m * pi * i * delta_x) / length));
		}
		exact[i] = T_sur + 2 * (T_init - T_sur) * series;
		v1.push_back(exact[i]);
	}

	/*cout << "The Exact Solution for this problem is: " << endl;
	cout << "--------------------------------------- " << endl;
	cout << "At Time t = " << fixed << setprecision(4) << output_time << " :" << endl;
	cout << "--------------------" << endl;
	for (int i = 0; i < nodes; i++) {
		cout << fixed << setprecision(4) << "Node #" << i << " at (x = " << (i * delta_x) << ") = " << exact[i] << endl;
	}
	cout << endl << endl;*/

	return v1;
}


Vector DuFort_Frankel_Explicit_Scheme(double** numerical, double delta_x, double delta_t, double output_time, int nodes, double r, double T_sur, double T_init) {
	Vector v2;

	for (int n = 0; n <= (output_time / delta_t); n++) {
		numerical[0][n] = T_sur;   // BC at 0 and all n (node #0)
		numerical[nodes - 1][n] = T_sur; // BC at 31 and all n (node #last_node)
	}

	for (int i = 1; i < (nodes - 1); i++) { // Initial Conditions #1
		numerical[i][0] = T_init; // IC at n = 0, for all i, except node node i = 0 and i = nodes
	}

	for (int i = 1; i < (nodes - 1); i++) { // Initial Conditions #2
		numerical[i][1] = r * numerical[i - 1][0] + (1.0 - (2.0 * r)) * numerical[i][0] + r * numerical[i + 1][0]; // IC at n = 1, for all i, except node node i = 0 and i = nodes
	}

	for (int n = 1; n < (output_time / delta_t); n++) {
		for (int i = 1; i < (nodes - 1); i++) {
			numerical[i][n + 1] = ((2.0 * r) / (1.0 + 2.0 * r)) * numerical[i - 1][n] + ((2.0 * r) / (1.0 + 2.0 * r)) * numerical[i + 1][n] + ((1.0 - 2.0 * r) / (1.0 + 2.0 * r)) * numerical[i][n - 1];
			
		}
	}

	int n = output_time / delta_t - 1;
	for (int i = 0; i < nodes; i++) {
		v2.push_back(numerical[i][n]);
	}

	/*cout << "The DuFort-Frankel Explicit Scheme Solution for this problem is: " << endl;
	cout << "---------------------------------------------------------------- " << endl;
	cout << "At Time t = " << output_time << " :" << endl;
	cout << "--------------------" << endl;
	int n = output_time / delta_t;
	for (int i = 0; i < nodes; i++) {
		cout << fixed << setprecision(4) << "Node #" << i << " at (x = " << (i * delta_x) << ") = " << numerical[i][n] << endl;
	}
	cout << endl << endl;*/

	return v2;
}


Vector Richardson_Explicit_Scheme(double** numerical, double delta_x, double delta_t, double output_time, int nodes, double r, double T_sur, double T_init) {
	Vector v3;

	for (int n = 0; n <= (output_time / delta_t); n++) {
		numerical[0][n] = T_sur;   // BC at 0 and all n (node #0)
		numerical[nodes - 1][n] = T_sur; // BC at 31 and all n (node #100)
	}

	for (int i = 1; i < (nodes - 1); i++) { // Initial Conditions #1
		numerical[i][0] = T_init; // IC at n = 0, for all i, except node node i = 0 and i = nodes
	}

	for (int i = 1; i < (nodes - 1); i++) { // Initial Conditions #2
		numerical[i][1] = r * numerical[i - 1][0] + (1.0 - (2.0 * r)) * numerical[i][0] + r * numerical[i + 1][0]; // IC at n = 1, for all i, except node node i = 0 and i = nodes
	}

	for (int n = 1; n < (output_time / delta_t); n++) {
		for (int i = 1; i < (nodes - 1); i++) {
			numerical[i][n + 1] = 2 * r * numerical[i - 1][n] - 4 * r * numerical[i][n] + 2 * r * numerical[i + 1][n] + numerical[i][n - 1];

		}
	}

	int n = output_time / delta_t - 1;
	for (int i = 0; i < nodes; i++) {
		v3.push_back(numerical[i][n]);
	}

	/*cout << "The Richardson Explicit Scheme Solution for this problem is: " << endl;
	cout << "------------------------------------------------------------ " << endl;
	cout << "At Time t = " << output_time << " :" << endl;
	cout << "--------------------" << endl;
	int n = output_time / delta_t;
	for (int i = 0; i < nodes; i++) {
		cout << fixed << setprecision(4) << "Node #" << i << " at (x = " << (i * delta_x) << ") = " << numerical[i][n] << endl;
	}
	cout << endl << endl;*/

	return v3;
}


Vector Laasonen_Simple_Implicit_Scheme(double** numerical, double delta_x, double delta_t, double output_time, int nodes, double r, double T_sur, double T_init) {
	Vector v4;

	for (int n = 0; n <= (output_time / delta_t); n++) {
		numerical[0][n] = T_sur;
		numerical[nodes - 1][n] = T_sur;
	}

	for (int i = 1; i < (nodes - 1); i++) {
		numerical[i][0] = T_init;					
	}
	
	double* lower_diag;
	double* main_diag;
	double* upper_diag;
	double* b;

	lower_diag = new double [nodes - 2];
	main_diag = new double [nodes - 2];
	upper_diag = new double [nodes - 2];
	b = new double [nodes - 2];

	for (int k = 0; k < nodes - 2; k++) {
		lower_diag[k] = 0;
		main_diag[k] = 0;
		upper_diag[k] = 0;
		b[k] = 0;
	}

	for (int k = 1; k < nodes - 2; k++) {
		lower_diag[k] = -r;
	}

	for (int k = 0; k < nodes - 2; k++) {
		main_diag[k] = 1 + 2 * r;
	}

	double* numerical_old;
	numerical_old = new double [nodes - 2];

	for (int n = 0; n < output_time / delta_t; n++) {
		for (int k = 0; k < nodes - 3; k++) {
			upper_diag[k] = -r;
		}

		for (int k = 0; k < nodes - 2; k++) {
			if (k == 0) {
				b[k] = numerical[k + 1][n] + r * T_sur;
			}
			else if (k == nodes - 3) {
				b[k] = numerical[k + 1][n] + r * T_sur;
			}
			else {
				b[k] = numerical[k + 1][n];
			}
		}

		TDMA_Solver(lower_diag, main_diag, upper_diag, b, nodes);
		for (int i = 0; i < nodes - 2; i++) {
			numerical[i + 1][n + 1] = b[i];
		}
	}
	
	int n = output_time / delta_t - 1;
	for (int i = 0; i < nodes; i++) {
		v4.push_back(numerical[i][n]);
	}

	/*cout << "The Laasonen Simple Implicit Scheme Solution for this problem is: " << endl;
	cout << "----------------------------------------------------------------- " << endl;
	cout << "At Time t = " << output_time << " :" << endl;
	cout << "--------------------" << endl;
	int n = output_time / delta_t;
	for (int i = 0; i < nodes; i++) {
			cout << fixed << setprecision(4) << "Node #" << i << " at (x = " << (i * delta_x) << ") = " << numerical[i][n] << endl;
	}
	cout << endl << endl;*/

	return v4;
}


Vector Crank_Nicholson_Implicit_Scheme(double** numerical, double delta_x, double delta_t, double output_time, int nodes, double r, double T_sur, double T_init) {
	Vector v5;

	for (int n = 0; n <= (output_time / delta_t); n++) {
		numerical[0][n] = T_sur;
		numerical[nodes - 1][n] = T_sur;
	}

	for (int i = 1; i < (nodes - 1); i++) {
		numerical[i][0] = T_init;
	}

	double* lower_diag;
	double* main_diag;
	double* upper_diag;
	double* b;

	lower_diag = new double[nodes - 2];
	main_diag = new double[nodes - 2];
	upper_diag = new double[nodes - 2];
	b = new double[nodes - 2];

	for (int k = 0; k < nodes - 2; k++) {
		lower_diag[k] = 0;
		main_diag[k] = 0;
		upper_diag[k] = 0;
		b[k] = 0;
	}
	
	for (int k = 1; k < nodes - 2; k++) {
		lower_diag[k] = -(r / 2);
	}

	for (int k = 0; k < nodes - 2; k++) {
		main_diag[k] = 1 + r;
	}
	
	double d = r / 2;
	double e = 1 - r;
	double f = r / 2;

	double* numerical_old;
	numerical_old = new double[nodes - 2];

	for (int n = 0; n < output_time / delta_t; n++) {
		for (int k = 0; k < nodes - 3; k++) {
			upper_diag[k] = -(r / 2);
		}

		for (int k = 0; k < nodes - 2; k++) {
			if (k == 0) {
				b[k] = d * numerical[k][n] + e * numerical[k + 1][n] + f * numerical[k + 2][n] + (r / 2) * T_sur;
			}
			else if (k == nodes - 3) {
				b[k] = d * numerical[k][n] + e * numerical[k + 1][n] + f * numerical[k + 2][n] + (r / 2) * T_sur;
			}
			else {
				b[k] = d * numerical[k][n] + e * numerical[k + 1][n] + f * numerical[k + 2][n];
			}
		}

		TDMA_Solver(lower_diag, main_diag, upper_diag, b, nodes);
		for (int i = 0; i < nodes - 2; i++) {
			numerical[i + 1][n + 1] = b[i];
		}
	}

	int n = output_time / delta_t - 1;
	for (int i = 0; i < nodes; i++) {
		v5.push_back(numerical[i][n]);
	}

	/*cout << "The Crank Nicholson Implicit Scheme Solution for this problem is: " << endl;
	cout << "----------------------------------------------------------------- " << endl;
	cout << "At Time t = " << output_time << " :" << endl;
	cout << "--------------------" << endl;
	int n = output_time / delta_t;
	for (int i = 0; i < nodes; i++) {
		cout << fixed << setprecision(4) << "Node #" << i << " at (x = " << (i * delta_x) << ") = " << numerical[i][n] << endl;
	}
	cout << endl << endl;*/

	return v5;
}


int TDMA_Solver(double* lower_diag, double* main_diag, double* upper_diag, double* b, int nodes) {
	nodes = nodes - 2;
	nodes--; 
	upper_diag[0] /= main_diag[0];
	b[0] /= main_diag[0];

	for (int i = 1; i < nodes; i++) {
		upper_diag[i] /= main_diag[i] - lower_diag[i] * upper_diag[i - 1];
		b[i] = (b[i] - lower_diag[i] * b[i - 1]) / (main_diag[i] - lower_diag[i] * upper_diag[i - 1]);
	}

	b[nodes] = (b[nodes] - lower_diag[nodes] * b[nodes - 1]) / (main_diag[nodes] - lower_diag[nodes] * upper_diag[nodes - 1]);

	for (int i = nodes; i-- > 0;) {
		b[i] -= upper_diag[i] * b[i + 1];
	}
	return 0;
}


void Print_Solution (Vector v1, Vector v2, Vector v3, Vector v4, Vector v5, int nodes, double delta_x) {
	ofstream outfile ("1D_Heat_Equation_Solution.csv");
	if (outfile.is_open()) {
		outfile << "x (m)" << "," << "Exact_Solution" << "," << "DuFort_Frankel_Explicit_Scheme" << "," << "Richardson_Explicit_Scheme" << "," << "Laasonen_Simple_Implicit_Scheme" << "," << "Crank_Nicholson_Implicit_Scheme" << endl;	
		for (int i = 0; i < nodes; i++) {
			outfile << fixed << setprecision(4) << (i * delta_x) << "," << v1[i] << "," << v2[i] << "," << v3[i] << "," << v4[i] << "," << v5[i] << endl;
		}
		outfile.close();
	}
}