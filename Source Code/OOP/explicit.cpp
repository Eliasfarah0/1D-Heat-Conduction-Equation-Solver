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


#include "explicit.h"
#include "implicit.h" // We use an Implicit object
#include <cmath>


// Default constructor
Explicit::Explicit() { // Will be initialised by the Analysis class
	deltat = 0;
	deltax = 0;
	D_value = 0;
	spaceDomain = 0;
	timeDomain = 0;
	t_surf = 0;
	t_init = 0;
}

// Get & set methods
void Explicit::setDeltat(double t) {
	deltat = t;
}
                                                                                           
void Explicit::setDeltax(double x) {
	deltax = x;
}

void Explicit::setD_value(double D) {
	D_value = D;
}

void Explicit::setSpaceDomain(int space) {
	spaceDomain = space;
}

void Explicit::setTimeDomain(int time) {
	timeDomain = time;
}

void Explicit::setT_surf(double Tsurf) {
	t_surf = Tsurf;
}

void Explicit::setT_init(double Tinit) {
	t_init = Tinit;
}

// Other methods
Vector Explicit::duFortSolve(int DufortFirstStepMethod) {
	double a = 2 * D_value * deltat / (deltax * deltax);
	int choice = DufortFirstStepMethod;
	Vector v1, v2, v3; // v1 == n - 1 // v2 == n // v3 == n + 1	

	switch (choice)
	{
		case 1: { // First Option: Use the FTCS scheme to get the solution at the first time step.
			// fill out the first vector i.e. initial temperature distribution along the space domain
			v1.push_back(t_surf); // node #0 (first node)
			for (int i = 1; i < spaceDomain; i++) { // nodes ranging from #1 to #619 (intermediate nodes)
				v1.push_back(t_init);
			}
			v1.push_back(t_surf); // node #620 (last node)

			// fill out the solution at the first time step
			v2.push_back(t_surf);
			for (int i = 1; i < spaceDomain; i++) {
				v2.push_back((a / 2) * v1[i - 1] + (1 - a) * v1[i] + (a / 2) * v1[i + 1]); // FTCS(forward in time, Central in space)
			}
			v2.push_back(t_surf);
			break;
		}
		case 2: { // Second Option: At t=0 every space node at 38C, and set the sides at 149C
			v1.push_back(t_init);
			for (int i = 1; i < spaceDomain; i++) {
				v1.push_back(t_init);
			}
			v1.push_back(t_init);

			v2.push_back(t_surf);
			for (int i = 1; i < spaceDomain; i++) {
				v2.push_back(t_init);
			}
			v2.push_back(t_surf);
			break;
		}
		case 3:	{ // Third Option: Use the laasonen simple implicit scheme for the first time step
			v1.push_back(t_surf);
			for (int i = 1; i < spaceDomain; i++) {
				v1.push_back(t_init);
			}
			v1.push_back(t_surf);

			// create an object Implicit
			Implicit laassonen;
			laassonen.setDeltat(deltat);
			laassonen.setDeltax(deltax);
			laassonen.setSpaceDomain(spaceDomain);
			laassonen.setTimeDomain(1);   // use this method only for the first time step
			laassonen.setD_value(D_value);
			laassonen.setT_init(t_init);
			laassonen.setT_surf(t_surf);
			v2 = laassonen.laasonenSolve();
			break;
		}
		case 4: { // Fourth Option: use FTCS but with a time step at 0.00001, so more likely stable than the first FTCS.
			v1.push_back(t_surf);
			for (int i = 1; i < spaceDomain; i++) {
				v1.push_back(t_init);
			}
			v1.push_back(t_surf);

			Vector v4;
			double b = 2 * D_value * 0.00001 / (deltax * deltax);
			for (int t = 0; t < deltat/0.00001; t++) { // adapt the number of iterration to stop at the first time step of Dufort-Frankel
				v4.push_back(t_surf);
				for (int i = 1; i < spaceDomain; i++) {
					v4.push_back((b / 2) * v1[i - 1] + (1 - b) * v1[i] + (b / 2) * v1[i + 1]);
				}
				v4.push_back(t_surf);
			}
			v2 = v4;
			v4.clear();
 			break;
		}
		default:
			std::cout << "ERROR! ENTER A VALUE OF 1, 2, 3 or 4 ONLY: ";
			break; 
	}

	// For the other time steps, use the classic DuFort-Frankel Scheme
	for (int t = 2; t < timeDomain; t++) {
		v3.push_back(t_surf);
		for (int i = 1; i < spaceDomain; i++) {
			v3.push_back(((1 - a) / (1 + a)) * v1[i] + (a / (1 + a)) * (v2[i + 1] + v2[i - 1]));
		}
		v3.push_back(t_surf);

		// stack management before the next loop
		v1 = v2;
		v2 = v3;
		v3 = {};
	}
	return v2; // Last Vector returned
}

Vector Explicit::richardsonSolve() {
	double a = 2 * D_value * deltat / (deltax * deltax);
	Vector v1, v2, v3;

	// fill out the initial vector
	v1.push_back(t_surf);
	for (int i = 1; i < spaceDomain; i++) {
		v1.push_back(t_init);
	}
	v1.push_back(t_surf);

	// use the FTCS method to get the solution at the first time step
	v2.push_back(t_surf);
	for (int i = 1; i < spaceDomain; i++) {
		v2.push_back((a / 2) * v1[i - 1] + (1 - a) * v1[i] + (a / 2) * v1[i + 1]);
	}
	v2.push_back(t_surf);

	// classic Richardson scheme to find out the other time step
	for (int t = 2; t < timeDomain; t++) {
		v3.push_back(t_surf);
		for (int i = 1; i < spaceDomain; i++) {
			v3.push_back(v1[i] + a * (v2[i + 1] - 2*v2[i] + v2[i - 1]));
		}
		v3.push_back(t_surf);

		// stack management before the next loop
		v1 = v2;
		v2 = v3;
		v3 = {};
	}
	return v2; // Last Vector returned
}