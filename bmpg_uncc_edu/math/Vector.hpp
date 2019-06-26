////////////////////////////////////////////////////////////////////////////////
//
// Copyright (c) 2008 by Mike Fairchild
// Author: Mike Fairchild, http://www.mikef.org, mjfairch@gmail.com
//
////////////////////////////////////////////////////////////////////////////////

#ifndef bmpg_uncc_edu_math_Vector_hpp
#define bmpg_uncc_edu_math_Vector_hpp

#include <cstddef>
#include <complex>

namespace bmpg_uncc_edu {
namespace math {

// Forward declaration for template friend operators.
// See: http://www.parashift.com/c++-faq-lite/templates.html#faq-35.16
template <typename T> class Vector;
template <typename T> std::ostream & operator<<(std::ostream & os, const Vector<T> & rhs);
template <typename T> std::istream & operator>>(std::istream & is, Vector<T> & rhs);

template <typename T>
class Vector
{
public:
	Vector();
	explicit Vector(size_t n,
			const T & value);
	
	explicit Vector(size_t n,
			T * data = NULL,
			bool pass_ownership = false);
	
	Vector(const Vector<T> & rhs);
	Vector<T> & operator=(const T & v);
	Vector<T> & operator=(const Vector<T> & rhs);
	~Vector();
	
	inline T * data() const { return data_; }
	T * release();
	void destroy();
	Vector<T> & assign(T * data,
			   size_t n,
			   bool pass_ownership = false);
	
	Vector<T> & reallocate(size_t n);
	Vector<T> & reallocate(size_t n,
			       const T & value);
	
	Vector<T> subVector(size_t i, size_t n); 

	inline size_t length() const { return length_; }
	T min() const;
	T max() const;
	T sum() const;
	inline T dot(const Vector<T> & v) const { return inner_product(v); }
	T inner_product(const Vector<T> & v) const;
	Vector<T> diff() const; // Returns Vector v such that v[i] = this[i+1]-this[i]
	T norm() const;
	T distance_to(const Vector<T> & v) const;
	T angle_to(const Vector<T> & v) const;

	Vector<T> & fill(T value);
	Vector<T> & normalize(T length = static_cast<T>(1));
	
	// Element access operators
	T element(size_t i) const;
	T & element(size_t i);
	T operator()(size_t i) const;
	T & operator()(size_t i);
	// Element-wise comparison operators
	bool operator==(const Vector<T> & v) const;
	inline bool operator!=(const Vector<T> & v) const { return !(operator==(v)); }
	// Element-wise Vector-scalar operations
	Vector<T> & operator+=(const T & v);
	Vector<T> & operator-=(const T & v);
	Vector<T> & operator*=(const T & v);
	Vector<T> & operator/=(const T & v);
	Vector<T> operator+(const T & v) const;
	Vector<T> operator-(const T & v) const;
	Vector<T> operator*(const T & v) const;
	Vector<T> operator/(const T & v) const;
	// Element-wise Vector-Vector operations
	Vector<T> & operator+=(const Vector<T> & v);
	Vector<T> & operator-=(const Vector<T> & v);
	Vector<T> & operator*=(const Vector<T> & v);
	Vector<T> & operator/=(const Vector<T> & v);
	Vector<T> operator+(const Vector<T> & v) const;
	Vector<T> operator-(const Vector<T> & v) const;
	Vector<T> operator*(const Vector<T> & v) const;
	Vector<T> operator/(const Vector<T> & v) const;
	
	friend std::ostream & operator<< <> (std::ostream & os, const Vector<T> & rhs);
	friend std::istream & operator>> <> (std::istream & os, Vector<T> & rhs);

protected:
	T * data_;
	bool own_data_;
	size_t length_;
};

} // namespace bmpg_uncc_edu::math
} // namespace bmpg_uncc_edu

#endif

