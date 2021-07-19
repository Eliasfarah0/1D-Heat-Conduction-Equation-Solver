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


#ifndef ANALYSIS_H
#define ANALYSIS_H
#include "explicit.h"  // we use Explicit objects in Analysis code
#include "implicit.h"  // we use Implicit objects in Analysis code


class Analysis {
	// Atributes
	private:
		double D_value, deltax, deltat, thickness, outputTime, t_surf, t_init; // respectively: diffusion coefficient, space step, time step, time which ,temperature of the sides, initial temperature 
		int DufortFirstStepMethod; // this integer will define witch method to use for getting the solution at the first time step of the Dufort-Frankel scheme
	
	public:
		// Default contructor
		Analysis(double Diff, double dx, double dt, double thick, double output_Time, double T_surf, double T_init, int choice_duFort);
		
		// Get & Set methods
		void setD_value(double D);
		void setDeltax(double x);
		void setDeltat(double t);
		void setThickness(double L);
		void setOutputTime(double outTime);
		void setT_surf(double Tsurf);
		void setT_init(double Tinit);
		void setDufortFirstStepMethod(int choice);
		
		// Methods
		Implicit initialiseImplicit(Implicit impl); // initialise the implicit object, especially define discret time and space domain
		Explicit initialiseExplicit(Explicit expl); // initialise the Explicit object, especially define discret time and space domain
		
		// write in a .csv file, the numerical solution at each node, using the Dufort-Frankel scheme, for a CONSTANT t chosen
		void printExplicit_duFordFrankel();
		
		// write in a .csv file, the numerical solution at each node, using the Richardson scheme, for a CONSTANT t chosen		
		void printExplicit_richardson();
		
		// write in a .csv file, the numerical solution at each node, using the Laasonen simple implicit scheme, for a CONSTANT t chosen
		void printImplicit_laasonen();
		
		// write in a .csv file, the numerical solution at each node, using the CrankNicolson scheme, for a CONSTANT t chosen
		void printImplicit_crankNicolson();
		
		// analytic solution for each nodes at the time we're looking for, in order to compare it with numerical solutions(errors)	
		Vector exact_solution();
		void print_exact_solution();
		
		// show the errors i.e. absolute difference betwteen numerical values and analytic values at each nodes
		Vector printErrors(int numerical_scheme);
	
		// write in a .csv file, the numerical solution at each time step until the duration chosen is reached, using a chosen numerical scheme, for a CONSTANT x
		void printTimeFunction(double positionToSee, double timeToSee, int numerical_scheme);
};
#endif
