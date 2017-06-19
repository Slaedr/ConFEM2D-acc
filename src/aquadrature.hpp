/** @file aquadrature.hpp
 * @brief Provides Gauss quadrature rules
 * @author Aditya Kashi
 * @date 2017-06-19
 */

#ifndef __AQUADRATURE_H
#define __AQUADRATURE_H

#ifndef __ADEFS_H
#include "adefs.hpp"
#endif

namespace afem2d {

/// Base class for defining quadrature rules
class QuadratureRule
{
protected:
	dShape shape;
	dushort ngauss;										///< Number of quadrature points
	dushort nPoly;										///< Degree of polynomial to integrate exactly
	dVector<a_real> gweights;
	dMatrix<a_real> ggpoints;
public:
	virtual void initialize(const int n_poly) = 0;
	virtual ~QuadratureRule() { }
	
	const dVector<a_real>& weights() const {
		return gweights;
	}

	const dMatrix<a_real>& points() const {
		return ggpoints;
	}

	dushort numGauss() const {
		return ngauss;
	}

	dShape getShape() const {
		return shape;
	}
};

/// 1D Gauss-Legendre quadrature
class Quadrature1D : public QuadratureRule
{
public:
	void initialize(const int n_poly);
};

class Quadrature2D : public QuadratureRule
{
public:
	virtual void initialize(const int n_poly) = 0;
};

/// Integration over the reference square
/** Note that currently, this is restricted to having the same number of quadrature points in the x- and y-directions.
 */
class Quadrature2DSquare : public Quadrature2D
{
public:
	void initialize(const int n_poly);
};

/// Integration over the reference triangle [(0,0), (1,0), (0,1)]
class Quadrature2DTriangle : public Quadrature2D
{
public:
	void initialize(const int n_poly);
};

} // end namespace acfd
#endif
