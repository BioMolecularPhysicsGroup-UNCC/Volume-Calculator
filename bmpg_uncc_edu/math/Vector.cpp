////////////////////////////////////////////////////////////////////////////////
//
// Copyright (c) 2008 by Mike Fairchild
// Author: Mike Fairchild, http://www.mikef.org, mjfairch@gmail.com
//
////////////////////////////////////////////////////////////////////////////////

#ifndef bmpg_uncc_edu_math_Vector_cpp
#define bmpg_uncc_edu_math_Vector_cpp

#include <complex>
#include <bmpg_uncc_edu/math/Vector.hpp>
#include <bmpg_uncc_edu/util/Exception.hpp>

namespace bmpg_uncc_edu {
namespace math {
using namespace std;

template <typename T>
Vector<T>::Vector() :
		  data_(NULL),
		  own_data_(false),
		  length_(0)
{
}

template <typename T>
Vector<T>::Vector(size_t n,
		  const T & value) : 
		  data_(NULL),
		  own_data_(false),
		  length_(0)
{
	reallocate(n, value);
}

template <typename T>
Vector<T>::Vector(size_t n,
		  T * d,
		  bool pass_ownership) : 
		  data_(NULL),
		  own_data_(false),
		  length_(0)
{
	if (d != NULL)
		assign(d, n, pass_ownership);
	else
		reallocate(n);
}

template <typename T>
Vector<T>::Vector(const Vector<T> & rhs) : 
	data_(NULL), own_data_(false), length_(0)
{
	operator=(rhs);
}

template <typename T>
Vector<T> & Vector<T>::operator=(const Vector<T> & rhs)
{
	if (this == &rhs)
		return *this;

	// attempt to allocate new memory
	reallocate(rhs.length_);
	register size_t i;
	for (i = 0; i < length_; ++i)
		data_[i] = rhs.data_[i];

	return *this;
}

template <typename T>
Vector<T>::~Vector()
{
	destroy();
}

template <typename T>
T * Vector<T>::release()
{
	if (data_ == NULL)
		throw bmpg_uncc_edu::util::Exception("Vector<T>::release(): data is NULL.", __FILE__, __LINE__);
	own_data_ = false;
	return data_;
}

template <typename T>
void Vector<T>::destroy()
{
	if (own_data_ && data_ != NULL)
		delete[] data_;
	data_ = NULL;
	length_ = 0;
	own_data_ = false;
}

template <typename T>
Vector<T> & Vector<T>::assign(T * d,
			      size_t n,
			      bool pass_ownership)
{
	destroy();
	data_ = d;
	length_ = n;
	own_data_ = pass_ownership;
	return *this;
}

template <typename T>
Vector<T> & Vector<T>::reallocate(size_t n)
{
	destroy();
	data_ = new T[n];
	if (data_ == NULL)
		throw bmpg_uncc_edu::util::Exception("Vector::Vector(n): Unable to allocate memory.", __FILE__, __LINE__);
	length_ = n;
	own_data_ = true;	//set own data true for new allocated memory  (Hui Wang)
	
	return *this;
}

template <typename T>
Vector<T> & Vector<T>::reallocate(size_t n,
				  const T & value)
{
	data_ = new T[n];
	if (data_ == NULL)
		throw bmpg_uncc_edu::util::Exception("Vector::Vector(n, T &): Unable to allocate memory.", __FILE__, __LINE__);
	length_ = n;
	fill(value);
	own_data_ = true;	//set own data true for new allocated memory (Hui Wang)
	
	return *this;
}

template <typename T>
Vector<T> Vector<T>::subVector(size_t i,
			       size_t n)
{
	if (i+n > length_)
		throw bmpg_uncc_edu::util::Exception("Vector::subVector(i,n): invalid arguments.");
	Vector<T> result(n);
	size_t j = i;
	for (; j-i < n; ++j)
		result.data_[j-i] = data_[j];
	return result;
}

template <typename T>
T Vector<T>::min() const
{
	if (length_ < 1)
		throw bmpg_uncc_edu::util::Exception("Vector::min(): Vector is empty.", __FILE__, __LINE__);
	T result = data_[0];
	register size_t i = 1;
	for (; i < length_; ++i) {
		if (data_[i] < result)
			result = data_[i];
	}
	return result;
}
template <> // Template specialization
std::complex<float> Vector<std::complex<float> >::min() const {
	throw bmpg_uncc_edu::util::Exception("Vector::min(): Undefined for complex Vectors.");
}
template <> // Template specialization
std::complex<double> Vector<std::complex<double> >::min() const {
	throw bmpg_uncc_edu::util::Exception("Vector::min(): Undefined for complex Vectors.");
}

template <typename T>
T Vector<T>::max() const
{
	if (length_ < 1)
		throw bmpg_uncc_edu::util::Exception("Vector::max(): Vector is empty.", __FILE__, __LINE__);
	T result = data_[0];
	register size_t i = 1;
	for (; i < length_; ++i) {
		if (data_[i] > result)
			result = data_[i];
	}
	return result;
}
template <> // Template specialization
std::complex<float> Vector<std::complex<float> >::max() const {
	throw bmpg_uncc_edu::util::Exception("Vector::max(): Undefined for complex Vectors.");
}
template <> // Template specialization
std::complex<double> Vector<std::complex<double> >::max() const {
	throw bmpg_uncc_edu::util::Exception("Vector::max(): Undefined for complex Vectors.");
}

template <typename T>
T Vector<T>::sum() const
{
	T result = static_cast<T>(0);
	if (length_ < 1)
		throw bmpg_uncc_edu::util::Exception("Vector::sum(): Vector is empty.", __FILE__, __LINE__);
	register size_t i = 0;
	for (; i < length_; ++i)
		result += data_[i];
	return result;
}

template <typename T>
T Vector<T>::inner_product(const Vector<T> & v) const
{
	if (length_ != v.length_)
		throw bmpg_uncc_edu::util::Exception("Vector::inner_product(): Vector size mismatch.", __FILE__, __LINE__);
	T result = static_cast<T>(0.0);
	register size_t i = 0;
	for (; i < length_; ++i)
		result += data_[i]*v.data_[i];
	return result;
}

template <typename T>
Vector<T> Vector<T>::diff() const
{
	if (length_ < 2)
		throw bmpg_uncc_edu::util::Exception("Vector::diff(): Vector needs at least two elements.", __FILE__, __LINE__);
	Vector<T> result(length_-1);
	register size_t i = 1;
	for (; i < length_; ++i)
		result.data_[i-1] = data_[i] - data_[i-1];
	return result;
}

template <typename T>
inline T Vector<T>::norm() const
{
	return std::sqrt(inner_product(*this));
}

template <typename T>
inline T Vector<T>::distance_to(const Vector<T> & v) const
{
	return ((v - *this).norm());
}

template <typename T>
T Vector<T>::angle_to(const Vector<T> & v) const
{
	T m = norm();
	T vm = v.norm();
	if (m == static_cast<T>(0) || vm == static_cast<T>(0))
		throw bmpg_uncc_edu::util::Exception("Vector::angle_between(): At least one Vector is the zero Vector, so angle is undefined.", __FILE__, __LINE__);
	
	return std::acos(inner_product(v)/(m*vm));
	// Select the right branch of acos (we want 0 <= angle < pi)
	//double a = acos(inner_product(v)/(m*vm));
	//if (a >= bmpg_uncc_edu::math::MathConst::pi)
	//	a = bmpg_uncc_edu::math::MathConst::pi - a;
	//return a;
}
template <> // Template specialization
std::complex<float> Vector<std::complex<float> >::angle_to(const Vector<std::complex<float> > & v) const 
{
	throw bmpg_uncc_edu::util::Exception("Vector::angle_between(): undefined for complex Vectors.");
}
template <> // Template specialization
std::complex<double> Vector<std::complex<double> >::angle_to(const Vector<std::complex<double> > & v) const 
{
	throw bmpg_uncc_edu::util::Exception("Vector::angle_between(): undefined for complex Vectors.");
}

template <typename T>
Vector<T> & Vector<T>::fill(T value)
{
	register size_t i = 0;
	for (; i < length_; ++i)
		data_[i] = value;
	return *this;
}

template <typename T>
Vector<T> & Vector<T>::normalize(T t)
{
	T mag = static_cast<T>(norm());
	if (mag == static_cast<T>(0))
		throw bmpg_uncc_edu::util::Exception("Vector::normalize(): Cannot normalize the zero Vector.", __FILE__, __LINE__);
	T scale = t/mag;
	operator*=(scale);
	return *this;
}

template <typename T>
inline T Vector<T>::element(size_t i) const
{
#ifdef bmpg_uncc_edu_DEBUG
	if (i >= length_)
		throw bmpg_uncc_edu::util::Exception("Vector<T>::element(): Bad index - indices must range 0 <= i <= Vector.length()-1", __FILE__, __LINE__);
#endif
	return data_[i];
}

template <typename T>
inline T & Vector<T>::element(size_t i)
{
#ifdef bmpg_uncc_edu_DEBUG
	if (i >= length_)
		throw bmpg_uncc_edu::util::Exception("Vector<T>::element(): Bad index - indices must range 0 <= i <= Vector.length()-1", __FILE__, __LINE__);
#endif
	return data_[i];
}

template <typename T>
inline T & Vector<T>::operator()(size_t i)
{
#ifdef bmpg_uncc_edu_DEBUG
	if (i >= length_)
		throw bmpg_uncc_edu::util::Exception("Vector<T>::operator(): Bad index - indices must range 0 <= i <= Vector.length()-1", __FILE__, __LINE__);
#endif
	return data_[i];
}

template <typename T>
inline T Vector<T>::operator()(size_t i) const
{
#ifdef bmpg_uncc_edu_DEBUG
	if (i >= length_)
		throw bmpg_uncc_edu::util::Exception("Vector<T>::operator(): Bad index - indices must range 0 <= i <= Vector.length()-1", __FILE__, __LINE__);
#endif
	return data_[i];
}

template <typename T>
Vector<T> & Vector<T>::operator+=(const Vector<T> & v)
{
	if (v.length_ != length_)
		throw bmpg_uncc_edu::util::Exception("Vector<T>::operator+=(): Vector size mismatch.", __FILE__, __LINE__);
	register size_t i = 0;
	for (; i < length_; ++i)
		data_[i] += v.data_[i];
	return *this;
}

template <typename T>
Vector<T> & Vector<T>::operator-=(const Vector<T> & v)
{
	if (v.length_ != length_)
		throw bmpg_uncc_edu::util::Exception("Vector<T>::operator-=(): Vector size mismatch.", __FILE__, __LINE__);
	register size_t i = 0;
	for (; i < length_; ++i)
		data_[i] -= v.data_[i];
	return *this;
}

template <typename T>
Vector<T> & Vector<T>::operator*=(const Vector<T> & v)
{
	if (v.length_ != length_)
		throw bmpg_uncc_edu::util::Exception("Vector<T>::operator*=(): Vector size mismatch.", __FILE__, __LINE__);
	register size_t i = 0;
	for (; i < length_; ++i)
		data_[i] *= v.data_[i];
	return *this;
}


template <typename T>
Vector<T> & Vector<T>::operator/=(const Vector<T> & v)
{
	if (v.length_ != length_)
		throw bmpg_uncc_edu::util::Exception("Vector<T>::operator/=(): Vector size mismatch.", __FILE__, __LINE__);
	register size_t i = 0;
	for (; i < length_; ++i)
		data_[i] /= v.data_[i];
	return *this;
}


template <typename T>
Vector<T> Vector<T>::operator+(const Vector<T> & v) const
{
	Vector<T> result(*this);
	result += v;
	return result;
}

template <typename T>
Vector<T> Vector<T>::operator-(const Vector<T> & v) const
{
	Vector<T> result(*this);
	result -= v;
	return result;
}

template <typename T>
Vector<T> Vector<T>::operator*(const Vector<T> & v) const
{
	Vector<T> result(*this);
	result *= v;
	return result;
}


template <typename T>
Vector<T> Vector<T>::operator/(const Vector<T> & v) const
{
	Vector<T> result(*this);
	result /= v;
	return result;
}


template <typename T>
Vector<T> & Vector<T>::operator+=(register const T & v)
{
	register size_t i = 0;
	for (; i < length_; ++i)
		data_[i] += v;
	return *this;
}

template <typename T>
Vector<T> & Vector<T>::operator-=(register const T & v)
{
	register size_t i = 0;
	for (; i < length_; ++i)
		data_[i] -= v;
	return *this;
}

template <typename T>
Vector<T> & Vector<T>::operator*=(register const T & v)
{
	register size_t i = 0;
	for (; i < length_; ++i)
		data_[i] *= v;
	return *this;
}


template <typename T>
Vector<T> & Vector<T>::operator/=(register const T & v)
{
	register size_t i = 0;
	for (; i < length_; ++i)
		data_[i] /= v;
	return *this;
}


template <typename T>
Vector<T> Vector<T>::operator+(const T & v) const
{
	Vector<T> result(*this);
	result += v;
	return result;
}

template <typename T>
Vector<T> Vector<T>::operator-(const T & v) const
{
	Vector<T> result(*this);
	result -= v;
	return result;
}

template <typename T>
Vector<T> Vector<T>::operator*(const T & v) const
{
	Vector<T> result(*this);
	result *= v;
	return result;
}


template <typename T>
Vector<T> Vector<T>::operator/(const T & v) const
{
	Vector<T> result(*this);
	result /= v;
	return result;
}


template <typename T>
bool Vector<T>::operator==(const Vector<T> & v) const
{
	if (v.length_ != length_)
		return false;
	register size_t i = 0;
	for (; i < length_; i++) {
		if (data_[i] != v.data_[i])
			return false;
	}
	return true;
}

template <typename T>
std::ostream & operator<<(std::ostream & os, const Vector<T> & rhs)
{
	// FIXME - different behavior if ios::binary set
	size_t i, n = rhs.length_;
	os << rhs.length_ << ' ';
	for (i = 0; i < n; ++i)
		os << rhs.data_[i] << ' ';
	return os;
}

template <typename T>
std::istream & operator>>(std::istream & is, Vector<T> & rhs)
{
	// FIXME - different behavior if ios::binary set
	size_t i = 0, n = 0;

	is >> n;
	rhs.reallocate(n);
	while (is && i < n)
		is >> rhs.data_[i++];
	if (i < n)
		throw bmpg_uncc_edu::util::Exception("operator>>(istream, Vector): stream ended prematurely.", __FILE__, __LINE__);

	return is;
}

// Explicit template instantiations
template class Vector<float>;
template class Vector<double>;
template class Vector<std::complex<float> >;
template class Vector<std::complex<double> >;

// Explicit template instanations of friend functions
template std::ostream & operator<<(std::ostream & os, const Vector<float> & rhs);
template std::ostream & operator<<(std::ostream & os, const Vector<double> & rhs);
template std::ostream & operator<<(std::ostream & os, const Vector<std::complex<float> > & rhs);
template std::ostream & operator<<(std::ostream & os, const Vector<std::complex<double> > & rhs);

template std::istream & operator>>(std::istream & is, Vector<float> & rhs);
template std::istream & operator>>(std::istream & is, Vector<double> & rhs);
template std::istream & operator>>(std::istream & is, Vector<std::complex<float> > & rhs);
template std::istream & operator>>(std::istream & is, Vector<std::complex<double> > & rhs);
} // namespace bmpg_uncc_edu::math
} // namespace bmpg_uncc_edu

#endif

