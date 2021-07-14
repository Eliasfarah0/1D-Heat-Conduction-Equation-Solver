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


#ifndef IMPLICIT_H 
#define IMPLICIT_H
#include "vector.h"


class Implicit {
	// Attributes
	private:
		Vector A, B, C; // Triadiagonal coefficient--A: coef for T_{i-1}, B:coef for T_{i}, C: coef for T_{i+1}
		int spaceDomain, timeDomain; // number of nodes for the space and time grids
		double deltat, deltax, D_value, t_surf, t_init; // respectively: time step, space step, diffusion coefficient, temperature of the sides, initial temperature.

	public:
		// Default contructor
		Implicit();
		
		// Get & Set methods
		void setDeltat(double t);
		void setDeltax(double x);
		void setD_value(double D);
		void setSpaceDomain(int space);
		void setTimeDomain(int time);
		void setT_surf(double Tsurf);
		void setT_init(double Tinit);
		
		// Thomas algorithm resolution: takes a vector T^{n} and return the vector T^{n+1}
		Vector thomas_algorithm(Vector v);
		
		// duFortSolve uses the Laasonen scheme to solve the heat equation.
		Vector laasonenSolve();
		
		// Crank Nicolson uses the Laasonen scheme to solve the heat equation.
		Vector crankNicolsonSolve();
};
#endif