////////////////////////////////////////////////////////////////////////////////
//
// Copyright (c) 2008 by Mike Fairchild
// Author: Mike Fairchild, http://www.mikef.org, mjfairch@gmail.com
//
////////////////////////////////////////////////////////////////////////////////

#ifndef bmpg_uncc_edu_math_R3Vector_hpp
#define bmpg_uncc_edu_math_R3Vector_hpp

#include <bmpg_uncc_edu/math/Vector.hpp>

namespace bmpg_uncc_edu {
namespace math {

template <typename T>
class R3Vector : public bmpg_uncc_edu::math::Vector<T>
{
public:
	R3Vector();
	R3Vector(const R3Vector<T> & rhs);
	R3Vector<T> & operator=(const bmpg_uncc_edu::math::Vector<T> & rhs);
	R3Vector<T> & operator=(const R3Vector<T> & rhs);
	explicit R3Vector(const bmpg_uncc_edu::math::Vector<T> & rhs);
	R3Vector(T x_, T y_, T z_);
	// Default operator= and destructor
	
	// Return the cross product between this Vector and the given Vector
	R3Vector<T> cross(const R3Vector<T> & v) const;
	// Return the Cartesian x coordinate of this Vector
	inline T x() const { return this->element(0); }
	// Return the Cartesian y coordinate of this Vector
	inline T y() const { return this->element(1); }
	// Return the Cartesian/Cylindrical z coordinate of this Vector
	inline T z() const { return this->element(2); }
	
	inline T at(int i) const {return this->element(i);}
	inline T& at(int i) {return this->element(i);}
	
	// Return the Spherical r coordinate of this Vector
	T r() const;
	// Return the Spherical theta coordinate of this Vector
	T theta() const;
	// Return the Spherical/Cylindrical phi coordinate of this Vector (branch cuts handled appropriately)
	T phi() const;
	// Return the Cylindrical rho coordinate of this Vector
	T rho() const;
	// Set the value of this Vector from provided Cartesian x,y,z coordinates
	R3Vector<T> & set_cartesian(T x_, T y_, T z_);
	// Set the value of this Vector from Spherical r,theta,phi coordinates
	R3Vector<T> & set_spherical(T r_, T theta_, T phi_);
	// Set the value of this Vector from Cylindrical rho,phi,z coordinates
	R3Vector<T> & set_cylindrical(T rho_, T phi_, T z_);
};

} // namespace bmpg_uncc_edu::math
} // namespace bmpg_uncc_edu

#endif


