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


#ifndef EXPLICIT_H
#define EXPLICIT_H
#include "vector.h" // We use vector objects as a data storage  


class Explicit {
	// Atributes
	private: 
		int spaceDomain, timeDomain; // number of nodes for the space and time grids
		double deltat, deltax, D_value, t_surf, t_init;  // respectively: time step, space step, diffusion coefficient, temperature of the sides, initial temperature.

	public:
		// Default contructor
		Explicit();
		
		// Get & Set methods
		void setDeltat(double t);
		void setDeltax(double x);
		void setD_value(double D);
		void setSpaceDomain(int space);
		void setTimeDomain(int time);
		void setT_surf(double Tsurf);
		void setT_init(double Tinit);
		
		// other Methods
		Vector duFortSolve(int DufortFirstStepMethod); // duFortSolve use the duFort Frankel scheme to solve the heat equation, the integer in parameter indicates which approximation will be carry out for the solution at the first time step
		Vector richardsonSolve(); // Same that duFortSolve method, but with the richarson method.
};
#endif