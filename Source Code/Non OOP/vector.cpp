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

/* *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *\
 * This .cpp file was provdided by our Instructors at Cranfield University *
\* *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  */


#include "vector.h"
#include <cmath>

// CONSTRUCTORS
/*=
* Default constructor (empty vector)
*/
Vector::Vector() : std::vector<double>() {}


/*
* Alternate constructor - creates a vector of a given size
*/
Vector::Vector(int Num) : std::vector<double>()
{
	// set the size
	(*this).resize(Num);

	// initialise with zero
	std::size_t i;
	for (i = 0; i < size(); i++) (*this)[i] = 0.0;
}

/*
* Copy constructor
*/
Vector::Vector(const Vector& copy) : std::vector<double>()
{
	(*this).resize(copy.size());
    // copy the data members (if vector is empty then num==0)
	std::size_t i;
    for (i=0; i<copy.size(); i++) (*this)[i]=copy[i]; 
}

/*
* accessor method - get the size
*/
int Vector::getSize() const
{
	return size();
}

// OVERLOADED OPERATORS
/*
* Operator= - assignment
*/
Vector& Vector::operator=(const Vector& copy)
{
	(*this).resize(copy.size());
	std::size_t i;
    for (i=0; i<copy.size(); i++) (*this)[i] = copy[i]; 
    return *this;
}


// COMPARISON
/*
* Operator== - comparison
*/
bool Vector::operator==(const Vector& v) const
{
	if (size() != v.size()) throw std::invalid_argument("incompatible vector sizes\n");
	std::size_t i;
	for (i = 0; i<size(); i++) 	if (fabs((*this)[i] - v[i]) > 1.e-07) { return false; }
	return true;
}

// NORMS
/*
* 1 norm
*/
double Vector::one_norm()
{
	int i;
	double sum = 0;
	int n = size();

	for ( i=0; i < n; i++ )
	    sum += fabs((*this)[i]);

	return sum;
}

/*
* 2 norm
*/
double Vector::two_norm()
{
	int i;
	double sum = 0;
	int n = size();

	for ( i=0; i < n; i++)
	    sum += fabs((*this)[i]*(*this)[i]);

	return (sqrt(sum));
}

/*
* uniform (infinity) norm
*/
double Vector::uniform_norm()
{
	int i;
	double max = 0;
	int n = size();

	for ( i=0; i<n; i++ )
	    if (max < fabs((*this)[i])) max = fabs((*this)[i]);
	return max;
}


// INPUT AND OUTPUT
/*
* keyboard input , user friendly
*/
std::istream& operator>>(std::istream& is, Vector& v)
{
	if (!v.size()) {
		int n;

		std::cout << "input the size for the vector" << std::endl;
		is >> n;
		//check input sanity
		if(n < 0) throw std::invalid_argument("read error - negative vector size");

		// prepare the vector to hold n elements
		v = Vector(n);
	}
	// input the elements
	std::cout << "input "<< v.size() <<" vector elements" << std::endl;
	std::size_t i;
	for (i=0; i<v.size(); i++) is >> v[i];

    // return the stream object
    return is;
}

/*
* file input - raw data, compatible with file writing operator
*/
std::ifstream& operator>>(std::ifstream& ifs, Vector& v) 
{
    int n;

    // read size from the file
    ifs >> n;
    //check input sanity
    if(n < 0) throw std::invalid_argument("file read error - negative vector size");

    // prepare the vector to hold n elements
    v = Vector(n);

    // input the elements
    for (int i=0; i<n; i++) ifs >> v[i];

    // return the stream object
    return ifs;
}

/*
* screen output, user friendly
*/
std::ostream& operator<<(std::ostream& os, const Vector& v)
{
    if (v.size() > 0) {
		std::size_t i;
        for (i=0; i<v.size(); i++) os << v[i]  << " ";
        os << std::endl;
    }
    else
    {
        os << "Vector is empty." << std::endl;
    }
    return os;
}

/*
* file output - raw data, comaptible with file reading operator
*/
std::ofstream& operator<<(std::ofstream& ofs, const Vector& v)
{
    //put vector size in first line (even if it is zero)
    ofs << v.size() << std::endl;
    //put data in second line (if size==zero nothing will be put)
	std::size_t i;
    for (i=0; i<v.size(); i++) ofs << v[i]  <<  " ";
    ofs << std::endl;
    return ofs;
}