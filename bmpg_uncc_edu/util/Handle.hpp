/**
   Copyright (c) 2008 by Hui Wang
   @author Hui Wang
   @date Jan 20, 2009
*/

#ifndef bmpg_uncc_edu_util_Handle_hpp
#define bmpg_uncc_edu_util_Handle_hpp

#include<boost/shared_ptr.hpp>

namespace bmpg_uncc_edu {
namespace util {

/**
	Handle to C++ objects, using boost::shared_ptr.
*/
template<class T>
class Handle
{
public:
	Handle(const boost::shared_ptr<T>& h = boost::shared_ptr<T>());
	Handle(T* h) : h_(h) {}
	
	template<class U>
	Handle(const Handle<U>& h);
	
	
	const boost::shared_ptr<T>& operator->() const;
//	const boost::shared_ptr<T>& operator*() const;
	T& operator*() const;
	
	template<class U>
	bool operator==(const Handle<U>& other) const;
	
	template<class U>
	bool operator!=(const Handle<U>& other) const;
	
	template<class U>
	bool operator<(const Handle<U>& other) const;
	
	//Handle<T>& operator=(const Handle<T>& other){h_ = other.value(); return *this;}
	
	template<class U>
	Handle<T>& operator=(const Handle<U>& other){h_ = other.value(); return *this;}

	boost::shared_ptr<T> value() const {return h_;}
	
	bool empty() const;
private:
	boost::shared_ptr<T> h_;
};

template<class T>
Handle<T>::Handle(const boost::shared_ptr<T>& h) :
		  h_(h)
{
	
}

template<class T>
template<class U>
Handle<T>::Handle(const Handle<U>& h)
{
	operator=(h);
}

template<class T>
bool Handle<T>::empty() const
{
	return !h_;
}

template<class T>
const typename boost::shared_ptr<T>& Handle<T>::operator->() const
{
	return h_;	
}
	
template<class T>
T& Handle<T>::operator*() const
//const typename boost::shared_ptr<T>& Handle<T>::operator*() const
{
	return *h_;
}

template<class T>
template<class U>
bool Handle<T>::operator==(const Handle<U>& other) const
{
	return (h_ == other.h_);
}


template<class T>
template<class U>
bool Handle<T>::operator!=(const Handle<U>& other) const
{
	return !(h_ == other.h_);
}

template<class T>
template<class U>
bool Handle<T>::operator<(const Handle<U>& other) const
{
	return (*h_ < *(other.h_));
}

}	//bmpg_uncc_edu::util
}	//bmpg_uncc_edu
#endif				//bmpg_uncc_edu_util_Handle_hpp

