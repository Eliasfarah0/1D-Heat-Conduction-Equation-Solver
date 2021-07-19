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


#include "analysis.h"
#include <fstream> // To write into a .CSV file
#include <cmath>
#include <iomanip>
using namespace std;


// Default contructor
Analysis::Analysis(double Diff, double dx, double dt, double L, double output_Time, double T_surf, double T_init, int choice_duFort) { // this->x == (*this).x
	(*this).D_value               = Diff;          // or this->D_value               = Diff
	(*this).deltax                = dx;            // or this->deltax                = dx
	(*this).deltat                = dt;            // or this->deltat                = dt
	(*this).thickness             = L;             // or this->thickness             = L
	(*this).outputTime            = output_Time;   // or this->outputTime            = output_Time
	(*this).t_surf                = T_surf;        // or this->t_surf                = T_surf
	(*this).t_init                = T_init;        // or this->t_init                = T_init
	(*this).DufortFirstStepMethod = choice_duFort; // or this->DufortFirstStepMethod = choice_duFort // This value can be changed in the main function
}

// Get & Set methods
void Analysis::setD_value(double D) {
	D_value = D;
}

void Analysis::setDeltax(double x) {
	deltax = x;
}

void Analysis::setDeltat(double t) {
	deltat = t;
}

void Analysis::setThickness(double L) {
	thickness = L;
}

void Analysis::setOutputTime(double outTime) {
	outputTime = outTime;
}

void Analysis::setT_surf(double Tsurf) {
	t_surf = Tsurf;
}

void Analysis::setT_init(double Tinit) {
	t_init = Tinit;
}

void  Analysis::setDufortFirstStepMethod(int choice) {
	DufortFirstStepMethod = choice;
}

// Methods
/* similar to:
type Analysis::initialiseImplicit(type impl) {
	...;
	return impl;
}*/
Implicit Analysis::initialiseImplicit(Implicit impl) {
	impl.setDeltat(deltat);
	impl.setDeltax(deltax);
	impl.setD_value(D_value);
	impl.setSpaceDomain(int(thickness / deltax)); // space domain set as the integer value of the thickness of the wall divided by the space step.
	impl.setTimeDomain(int(outputTime / deltat)); // space domain set as the integer value of the thickness of the wall divided by the time step.
	impl.setT_surf(t_surf);
	impl.setT_init(t_init);
	return impl;
}

Explicit Analysis::initialiseExplicit(Explicit expl) {
	expl.setDeltat(deltat);
	expl.setDeltax(deltax);
	expl.setD_value(D_value);
	expl.setSpaceDomain(int(thickness / deltax)); // space domain set as the integer value of the thickness of the wall divided by the space step.
	expl.setTimeDomain(int(outputTime / deltat)); // space domain set as the integer value of the thickness of the wall divided by the time step.
	expl.setT_surf(t_surf);
	expl.setT_init(t_init);
	return expl;
}

void Analysis::printExplicit_duFordFrankel() {
	//creation and initialisation of an Explicit object
	Explicit solver;
	solver = initialiseExplicit(solver);
	/* solver = initialiseExplicit(solver); ---> this replace all the following arguments:
	Explicit solver;
	solver.setDeltat(deltat);
	solver.setDeltax(deltax);
	solver.setD_value(D_value);
	solver.setSpaceDomain(int(thickness / deltax));
	solver.setTimeDomain(int(outputTime / deltat));
	solver.setT_surf(t_surf);
	solver.setT_init(t_init);*/
	Vector v1 = solver.duFortSolve(DufortFirstStepMethod);
	ofstream outfile("duFortFrankel.csv");
	if (outfile.is_open()) {
		outfile << "x (m)" << "," << "T (K)" << endl;
		for (int i = 0; i < v1.size(); i++) {
			outfile << fixed << setprecision(4) << (i * deltax) << "," << v1[i] << endl;
		}
		outfile.close();
	}
}

void Analysis::printExplicit_richardson() {
	Explicit solver;
	solver = initialiseExplicit(solver);
	Vector v1 = solver.richardsonSolve();
	ofstream outfile("Richardson.csv");
	if (outfile.is_open()) {
		outfile << "x (m)" << "," << "T (K)" << endl;
		for (int i = 0; i < v1.size(); i++) {
			outfile << fixed << setprecision(4) << (i * deltax) << "," << v1[i] << endl;
		}
		outfile.close();
	}
}

void Analysis::printImplicit_laasonen() {
	Implicit solver;
	solver = initialiseImplicit(solver);
	Vector v1 = solver.laasonenSolve();
	ofstream outfile("laasonen.csv");
	if (outfile.is_open()) {
		outfile << "x (m)" << "," << "T (K)" << endl;
		for (int i = 0; i < v1.size(); i++) {
			outfile << fixed << setprecision(4) << (i * deltax) << "," << v1[i] << endl;
		}
		outfile.close();
	}
}

void Analysis::printImplicit_crankNicolson() {
	Implicit solver;
	solver = initialiseImplicit(solver);
	Vector v1 = solver.crankNicolsonSolve();
	ofstream outfile("crankNicolson.csv");
	if (outfile.is_open()) {
		outfile << "x (m)" << "," << "T (K)" << endl;
		for (int i = 0; i < v1.size(); i++) {
			outfile << fixed << setprecision(4) << (i * deltax) << "," << v1[i] << endl;
		}
		outfile.close();
	}
}

Vector Analysis::exact_solution() {
	double pi = 3.1415926535;
	Vector v1;
	int acc = 100;
	int spaceDomain = int(thickness / deltax);
	for (int i = 0; i < spaceDomain + 1; i++) {
		double sum = 0.0;
		for (int j = 1; j < acc; j++) {
			sum += exp(-D_value * pow((j * pi / (thickness)), 2) * outputTime) * ((1 - pow(-1, j)) / (j * pi)) * sin(j * pi * (i * deltax) / (thickness)); // Analytical (Exact) Solution 
		}
		v1.push_back(t_surf + (2 * (t_init - t_surf) * sum));
	}
	return v1;
}

void Analysis::print_exact_solution() {
	Vector v1 = exact_solution();
	///write the exact solution and return it as well.
	ofstream outfile("exact_solution.csv");
	if (outfile.is_open()) {
		outfile << "x (m)" << "," << "T (K)" << endl;
		for (int i = 0; i < v1.size(); i++) {
			outfile << fixed << setprecision(4) << (i * deltax) << "," << v1[i] << endl;
		}
		outfile.close();
	}
}

Vector Analysis::printErrors(int numerical_scheme) {
	Vector v1 = exact_solution();
	Vector v2;
	Vector errors(v1.size());
	string file;

	Explicit expl;
	expl = initialiseExplicit(expl);
	Implicit impl;
	impl = initialiseImplicit(impl);
	
	// for the numerical numerical_scheme chosen, write the errors in a .csv file, and return them as well
	switch (numerical_scheme) {
		case 1: {
			v2 = expl.duFortSolve(DufortFirstStepMethod);
			file = "errors_duFort.csv";
			break;
		}
		case 2: {
			v2 = expl.richardsonSolve();
			file = "errors_richardson.csv";
			break;
		}
		case 3: {
			v2 = impl.laasonenSolve();
			file = "errors_laasonen.csv";
			break;
		}
		case 4: {
			v2 = impl.crankNicolsonSolve();
			file = "errors_crankNicolson.csv";
			break;
		}
	}

	ofstream outfile(file);
	if (outfile.is_open()) {
		outfile << "x (m)" << "," << "Error" << endl;
		for (int i = 0; i < v1.size(); i++) {
			errors[i] = abs(v1[i] - v2[i]);
			outfile << (i * deltax) << "," << errors[i] << endl;
		}
		outfile.close(); 
	}
	return errors;
}

void Analysis::printTimeFunction(double positionToSee, double timeToSee, int numerical_scheme) {
	// for the numerical_scheme chosen (int numerical_scheme), evolution of the temperature at the node int(positionToSee/delta x) in time until t=timeToSee
	string file;
	int space = int(positionToSee / deltax);
	int time = int(timeToSee / deltat);
	Vector v1 = {};
	Explicit expl;
	expl = initialiseExplicit(expl);
	Implicit impl;
	impl = initialiseImplicit(impl);

	// check wether the node we are looking at is inside the domain
	if (positionToSee <= thickness) {
		for (int t = 0; t <= time; t++) {
			switch (numerical_scheme) {
				case 1: {
					file = "timeFunction_duFort.csv";
					expl.setTimeDomain(t);
					v1.push_back(expl.duFortSolve(DufortFirstStepMethod)[space]);
					break;
				}
				case 2: {
					file = "timeFunction_richardson.csv";
					expl.setTimeDomain(t);
					v1.push_back(expl.richardsonSolve()[space]);
					break;
				}
				case 3: {
					file = "timeFunction_laasonen.csv";
					impl.setTimeDomain(t);
					v1.push_back(impl.laasonenSolve()[space]);
					break;
				}
				case 4: {
					file = "timeFunction_crankNicolson.csv";
					impl.setTimeDomain(t);
					v1.push_back(impl.crankNicolsonSolve()[space]);
					break;
				}
			}
		}

		ofstream outfile(file);
		if (outfile.is_open()) {
			outfile << "At x = " << positionToSee << endl;
			outfile << "t (s)" << "," << "T (K)" << endl;
			for (int i = 0; i < v1.size(); i++) {
				outfile << (i * deltat) << fixed << setprecision(4) << "," << v1[i] << endl;
			}
			outfile.close(); 
		}
	}
	else
		cout << "The value of x chosen is out of borders!" << endl;
}
