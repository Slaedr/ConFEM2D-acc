/** @file adefs.hpp
 * @brief Global includes and definitions
 */

#ifndef __ADEFS_H
#define __ADEFS_H

#ifndef VIENNACL_SCALAR_HPP_
#include <viennacl/scalar.hpp>
#endif
#ifndef VIENNACL_VECTOR_HPP_
#include <viennacl/vector.hpp>
#endif
#ifndef VIENNACL_MATRIX_HPP_
#include <viennacl/matrix.hpp>
#endif

#define  PI   3.14159265358979323846
#define SQRT2 1.4142135623730950
#define SQRT3 1.73205080756887729353
#define ZERO_TOL 2.2e-16

#define LINE 0
#define TRIANGLE 1
#define QUADRANGLE 2

namespace afem2d 
{
	typedef double a_real;
	typedef unsigned int a_int;

	typedef unsigned short Shape;
	typedef unsigned short ushort;

	/// Device scalar encoding shape of element
	typedef viennacl::scalar<unsigned short> dShape;
	/// Device scalar for other small counting integers
	typedef viennacl::scalar<unsigned short> dushort;
	/// Device scalar for indices
	typedef viennacl::scalar<a_int> da_int;
	/// Device scalar for computations
	typedef viennacl::scalar<a_real> da_real;
	/// Device vector
	template<typename T>
	using dVector = viennacl::vector<T>;
	/// Device matrix
	template<typename T>
	using dMatrix = viennacl::matrix<T>;
}

#endif
