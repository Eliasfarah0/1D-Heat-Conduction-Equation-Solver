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
 *  This .h file was provdided by our Instructors at Cranfield University  *
\* *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  */


#ifndef VECTOR_H //Include guard
#define VECTOR_H

#include <iostream> // Generic IO operations
#include <fstream>  // File IO operations
#include <stdexcept> // provides exceptions
#include <vector>  // std vector upon which our Vector is based


/**
*  A vector class for data storage of a 1D array of doubles
*  \n The implementation is derived from the standard container vector std::vector
*  \n We use private inheritance to base our vector upon the library version whilst
*  \n allowing us to expose only those base class functions we wish to use - in this
*  \n case the array access operator []
*
* The Vector class provides:
* \n-basic constructors for creating vector obcjet from other vector object,
or by creating empty vector of a given size,
* \n-input and oput operation via >> and << operators using keyboard or file
* \n-basic operations like access via [] operator, assignment and comparision
*/
class Vector : public std::vector<double> {
	typedef std::vector<double> vec;
public:

	// CONSTRUCTORS
	/**
	* Default constructor.  Intialize an empty Vector object
	* @see Vector(int Num)
	* @see Vector(const Vector& v)
	*/
    Vector(); // default constructor


	/**
	* Explicit alterative constructor takes an intiger.
	* it is explicit since implicit type conversion int -> vector doesn't make sense
	* Intialize Vector object of size Num
	* @see Vector()
	* @see Vector(const Vector& v)
	* @exception invalid_argument ("vector size negative")
	*/
    explicit Vector(int Num /**< int. Size of a vector */);

	/**
	* Copy constructor takes an Vector object reference.
	* Intialize Vector object with another Vector object
	* @see Vector()
	* @see Vector(int Num)
	*/
    Vector(const Vector& v); 


	// OVERLOADED OPERATORS
	/**
	* Overloaded assignment operator
	* @see operator==(const Vector& v)const
	* @param v Vector to assign from
	* @return the object on the left of the assignment
	*/
    Vector& operator=(const Vector& v  /**< Vecto&. Vector to assign from */);

	/**
	* Overloaded comparison operator
	* returns true if vectors are the same within a tolerance (1.e-07)
	* @see operator=(const Vector& v)
	* @see operator[](int i)
	* @see operator[](int i)const
	* @return bool. true or false
	* @exception invalid_argument ("incompatible vector sizes\n")
	*/
	bool operator==(const Vector& v /**< Vector&. vector to compare */
		) const;


	// ACCESSOR METHODS
	/** Normal get method that returns integer, the size of the vector
	* @return int. the size of the vector
	*/
	int getSize() const;  


	// THREE NORMS
	/**
	* Normal public method that returns a double.
	* It returns L1 norm of vector
	* @see two_norm()const
	* @see uniform_norm()const
	* @return double. vectors L1  norm
	*/
	double one_norm();

	/**
	* Normal public method that returns a double.
	* It returns L2 norm of vector
	* @see one_norm()const
	* @see uniform_norm()const
	* @return double. vectors L2  norm
	*/
	double two_norm();

	/**
	* Normal public method that returns a double.
	* It returns L_max norm of vector
	* @see one_norm()const
	* @see two_norm()const
	* @exception out_of_range ("vector access error")
	* vector has zero size
	* @return double. vectors Lmax  norm
	*/
	double uniform_norm();


	// KEYBOARD/SCREEN INPUT AND OUTPUT
	/**
	* Overloaded istream >> operator. Keyboard input
	* if vector has size user will be asked to input only vector values
	* if vector was not initialized user can choose vector size and input it values
	* @see operator>>(std::ifstream& ifs, Vector& v)
	* @see operator<<(std::ostream& os, const Vector& v)
	* @see operator<<(std::ofstream& ofs, const Vector& v)
	* @return std::istream&. the input stream object is
	* @exception std::invalid_argument ("read error - negative vector size");
	*/
    friend std::istream& operator >> (std::istream& is, /**< keyboard input straem. For user input */
		Vector& v /**< Vector&. vector to write to */
		);


	/**
	* Overloaded ifstream << operator. Display output.
	* @see operator>>(std::istream& is, Vector& v)
	* @see operator>>(std::ifstream& ifs, Vector& v)
	* @see operator<<(std::ofstream& ofs, const Vector& v)
	* @return std::ostream&. the output stream object os
	*/
	friend std::ostream& operator<<(std::ostream& os, /**< output file stream */
		const Vector& v /**< vector to read from */
		);

	/**
	* Overloaded ifstream >> operator. File input
	* the file output operator is compatible with file input operator,
	* ie. everything written can be read later.
	* @see operator>>(std::istream& is, Vector& v)
	* @see operator<<(std::ostream& os, const Vector& v)
	* @see operator<<(std::ofstream& ofs, const Vector& v)
	* @return ifstream&. the input ifstream object ifs
	* @exception std::invalid_argument ("file read error - negative vector size");
	*/
	friend std::ifstream& operator >> (std::ifstream& ifs, /**< input file straem. With opened matrix file */
		Vector& v /**< Vector&. vector to write to */
		);


	/**
	* Overloaded ofstream << operator. File output.
	* the file output operator is compatible with file input operator,
	* ie. everything written can be read later.
	* @see operator>>(std::istream& is, Vector& v)
	* @see operator>>(std::ifstream& ifs, Vector& v)
	* @see operator<<(std::ostream& os, const Vector& v)
	* @return std::ofstream&. the output ofstream object ofs
	*/
	friend std::ofstream& operator<<(std::ofstream& ofs, /**< outputfile stream. With opened file */
		const Vector& v /**< Vector&. vector to read from */
		);
};

#endif