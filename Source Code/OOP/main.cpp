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
using namespace std;


int main()
{
	// Problem Variables (Diffusivity, DeltaX, DeltaT, Thickness, OutputTime, Tsurf, Tinit, duFortFirstStepMethod);
	Analysis HeatEquation(93, 0.05, 0.01, 31, 0.5, 149, 38, 1);

	// Exact Analytical Solution for the 1D Heat Equation
	HeatEquation.print_exact_solution();

	// Numerical Solution Using 4 Different Numerical Schemes
	HeatEquation.printExplicit_duFordFrankel();
	HeatEquation.printExplicit_richardson();
	HeatEquation.printImplicit_laasonen();
	HeatEquation.printImplicit_crankNicolson();
	
	// Absolute Errors for All Numerical Schemes With Respect to the Exact Analytical Solution
	// HeatEquation.printErrors(numerical_scheme);
	HeatEquation.printErrors(1);
	HeatEquation.printErrors(2);
	HeatEquation.printErrors(3);
	HeatEquation.printErrors(4);
	
	// Numerical solution at a specified location, until a specified time and using a specified numerical scheme
	// HeatEquation.printTimeFunction(positionToSee, timeToSee, numerical_scheme);
	HeatEquation.printTimeFunction(5, 0.5, 1);
	HeatEquation.printTimeFunction(10, 0.5, 2);
	HeatEquation.printTimeFunction(15.5, 0.5, 3);
	HeatEquation.printTimeFunction(20, 0.5, 4);

	cout << "Computation Completed!" << endl;
	system("PAUSE");
	return 0;
}
