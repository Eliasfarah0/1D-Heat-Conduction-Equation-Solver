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


#include "implicit.h"
#include <cmath>


// Default constructor
Implicit::Implicit() {
	deltat = 0;
	deltax = 0;
	D_value = 0;
	spaceDomain = 0;
	timeDomain = 0;
	t_surf = 0;
	t_init = 0;
	A = {};
	B = {};
	C = {};
}

// Get & set methods
void Implicit::setDeltat(double t) {
	deltat = t;
}

void Implicit::setDeltax(double x) {
	deltax = x;
}

void Implicit::setD_value(double D) {
	D_value = D;
}

void Implicit::setSpaceDomain(int space) {
	spaceDomain = space;
}

void Implicit::setTimeDomain(int time) {
	timeDomain = time;
}

void Implicit::setT_surf(double Tsurf) {
	t_surf = Tsurf;
}

void Implicit::setT_init(double Tinit) {
	t_init = Tinit;
}

Vector Implicit::thomas_algorithm(Vector d) {
	// fill out the three diagonals
	Vector a = A; // lower_diagonal
	Vector b = B; //  main_diagonal
	Vector c = C; // upper_diagonal

	int n = A.size();

	// row 0 divided by b
	c[0] /= b[0];
	d[0] /= b[0];
	b[0] = 1;

	for (int i = 1; i < n; i++) {// For all the rows left
		// row i -= row_{i-1} * a_{i}
		b[i] -= c[i - 1] * a[i]; 
		d[i] -= d[i - 1] * a[i];
		a[i] = 0; // a_{i}-1*a_{i}=0

		// row_{i} divided by b_{i}
		c[i] /= b[i];
		d[i] /= b[i];
		b[i] = 1; 
	}
	// d_{n} is stored as a solution.
	// all the d_{i} coefficients are set thank to the simplified sysem of equations
	for (int i = n - 2; i >= 0; i--) {
		d[i] -= (c[i] * d[i + 1]);
	}
	return d;
}

Vector Implicit::laasonenSolve() {
	Vector D;
	double a = D_value * (deltat / (deltax * deltax));

	// Fill out the three diagonals
	// Fill out the initial vector
	A.push_back(0); // lower_diagonal
	B.push_back(1); //  main_diagonal
	C.push_back(0); // upper_diagonal
	D.push_back(t_surf); // Solution_diagonal
	for (int i = 1; i < spaceDomain; i++) {
		A.push_back(-a);
		B.push_back(1 + (2 * a));
		C.push_back(-a);
		D.push_back(t_init);
	}
	A.push_back(0);
	B.push_back(1);
	C.push_back(0);   
	D.push_back(t_surf);

	// resset of D at each tme step
	for (int t = 1; t < timeDomain; t++) {
		D = thomas_algorithm(D);
	}

	// clear the diagonal in case of an other call of this method without initialisation
	A.clear();
	B.clear();
	C.clear();
	return D;
}

Vector Implicit::crankNicolsonSolve() {
	Vector D, init;
	double a = D_value * (deltat / (deltax * deltax));

	// Fill out the initial vector
	init.push_back(t_surf);
	for (int i = 1; i < spaceDomain; i++) {
		init.push_back(t_init);
	}
	init.push_back(t_surf);

	// size-2 for the diagonals
	A.push_back(0);
	B.push_back(a + 1);
	C.push_back(-a*0.5);
	for (int i = 1; i < spaceDomain-2; i++) {
		A.push_back(-a * 0.5);
		B.push_back(a + 1);
		C.push_back(-a * 0.5);
	}
	A.push_back(-a*0.5);
	B.push_back(1 + a);
	C.push_back(0);

	// reduce the size of the system N to N-2. Keep taking in count the boundary conditions by added to the right hand member of the system the values erased from the reduction, so -a*149 to d[1] and -c*149 to d[N-2]
	for(int t = 1; t < timeDomain; t++) {
		for (int i = 0; i < spaceDomain-1; i++) {
			if (i == 0)
				D.push_back((a / 2) * init[i] + (1 - a) * init[i + 1] + (a / 2) * init[i + 2] + (a / 2) * t_surf); 
			else if (i == spaceDomain - 2)
				D.push_back((a / 2) * init[i] + (1 - a) * init[i + 1] + (a / 2) * init[i + 2] + (a / 2) * t_surf);
			else 
				D.push_back((a / 2) * init[i] + (1 - a) * init[i + 1] + (a / 2) * init[i + 2]);
		}

		D = thomas_algorithm(D);

		// reset init by filling it with D (smaller size) and keep the boundarie conditions (init[0])
		for (int i = 0; i < D.size();i++) {
			init[i + 1] = D[i];
		}
		D.clear();
	}
	// clear the diagonal in case of an other call of this method without initialisation
	A.clear();
	B.clear();
	C.clear();
	return init;
}