////////////////////////////////////////////////////////////////////////////////
//
// Copyright (c) 2008 by Mike Fairchild
// Author: Mike Fairchild, http://www.mikef.org, mjfairch@gmail.com
//
////////////////////////////////////////////////////////////////////////////////

#ifndef bmpg_uncc_edu_math_R3Vector_cpp
#define bmpg_uncc_edu_math_R3Vector_cpp

#include <cmath>
#include <bmpg_uncc_edu/math/R3Vector.hpp>
#include <bmpg_uncc_edu/math/MathConstants.hpp>
#include <bmpg_uncc_edu/util/Exception.hpp>


#define drand48() (((double)rand())/((double)RAND_MAX)) //jenny

namespace bmpg_uncc_edu {
namespace math {
using namespace bmpg_uncc_edu::util;

template <typename T>
R3Vector<T>::R3Vector() : bmpg_uncc_edu::math::Vector<T>(3, static_cast<T>(0))
{
}

template <typename T>
R3Vector<T>::R3Vector(const R3Vector<T> & rhs) : bmpg_uncc_edu::math::Vector<T>(rhs.length(), static_cast<T>(0))
{
	operator=(rhs);
}

// This copy constructor is necessary and is NOT redundant.  For example, when 
// working with an object of type R3Vector and invoking a superclass method that
// returns a Vector<double> &, such as the method normalize(), for example, 
// we'd get compiler errors if it weren't for this copy constructor.
template <typename T>
R3Vector<T>::R3Vector(const bmpg_uncc_edu::math::Vector<T> & rhs) : bmpg_uncc_edu::math::Vector<T>(3, static_cast<T>(0))
{
	operator=(rhs);
}

template <typename T>
R3Vector<T> & R3Vector<T>::operator=(const bmpg_uncc_edu::math::Vector<T> & rhs)
{
	if (this == &rhs)
		return *this;
	
	if (rhs.length() != 3)
		throw Exception("R3Vector::R3Vector(R3Vector<T> &): Input Vector must have size() == 3.");
	
	set_cartesian(rhs(0), rhs(1), rhs(2));
	return *this; 
}

template <typename T>
R3Vector<T> & R3Vector<T>::operator=(const R3Vector<T> & rhs)
{
	if (this == &rhs)
		return *this;
	
	if (rhs.length() != 3)
		throw Exception("R3Vector::R3Vector(R3Vector<T> &): Input Vector must have size() == 3.");
	
	set_cartesian(rhs(0), rhs(1), rhs(2));
	return *this;
}

template <typename T>
R3Vector<T>::R3Vector(T x_, T y_, T z_) : bmpg_uncc_edu::math::Vector<T>(3, static_cast<T>(0))
{
	set_cartesian(x_, y_, z_); 
}

template <typename T>
R3Vector<T> R3Vector<T>::cross(const R3Vector<T> & v) const
{
	return R3Vector(y()*v.z()-v.y()*z(), v.x()*z()-x()*v.z(), x()*v.y()-v.x()*y());
}

template <typename T>
T R3Vector<T>::r() const
{
	return this->norm();
}

template <typename T>
T R3Vector<T>::theta() const
{
	R3Vector zhat(0,0,1);
//	return angle_to(zhat);        
	return this->angle_to(zhat); //jenny
}

template <typename T>
T R3Vector<T>::phi() const
{
	T x_ = x();
	T y_ = y();

	if (x_ == 0 && y_ == 0)
		throw "R3Vector::phi(): phi is undefined when x = y = 0";
	
	// phi is the angle between the projection of this Vector in the x-y plane
	// and the unit Vector in the x direction
	R3Vector xhat(1,0,0);
	R3Vector proj(x_, y_, 0);
	T ang = proj.angle_to(xhat);
	
	// We must select the right branch of acos depending on signs
	if (x_ < 0 && y_ < 0)
		return bmpg_uncc_edu::math::MathConstants::pi/2+ang;
	if (x_ > 0 && y_ < 0)
		return 2*bmpg_uncc_edu::math::MathConstants::pi-ang;
	if (x_ == 0 && y_ < 0)
		return bmpg_uncc_edu::math::MathConstants::pi+ang;
	return ang;
}

template <typename T>
T R3Vector<T>::rho() const
{
	return static_cast<T>(sqrt(pow(x(),2)+pow(y(),2)));
}

// NOTE: If you're wondering why I have all the this-> pointers, it's where I'm
// using something from the base class that doesn't depend on the template type.
// Since specializations of the base class might override these members differently,
// we must use a type-dependent way of accessing them, either throgh vec<T>::
// or this->, which is also type-dependent.
// For reference, see: http://gcc.gnu.org/onlinedocs/gcc/Name-lookup.html
template <typename T>
R3Vector<T> & R3Vector<T>::set_cartesian(T x_, T y_, T z_)
{
	this->element(0) = x_;
	this->element(1) = y_;
	this->element(2) = z_;
	return *this;
}

template <typename T>
R3Vector<T> & R3Vector<T>::set_spherical(T r_, T theta_, T phi_)
{
	set_cartesian(r_*sin(theta_)*cos(phi_), r_*sin(theta_)*sin(phi_), r_*cos(theta_));  
	return *this;
}

template <typename T>
R3Vector<T> & R3Vector<T>::set_cylindrical(T rho_, T phi_, T z_)
{
	set_cartesian(rho_*cos(phi_), rho_*sin(phi_), z_); 
	return *this;
}

// Explicit template instantiations
template class R3Vector<float>;
template class R3Vector<double>;

} // namespace bmpg_uncc_edu::math
} // namespace bmpg_uncc_edu

#endif


