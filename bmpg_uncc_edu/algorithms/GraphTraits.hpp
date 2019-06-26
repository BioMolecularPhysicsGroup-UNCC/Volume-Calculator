/**
   Copyright (c) 2008 by Hui Wang
   @author Hui Wang
   @date Jan 20, 2009
*/


#ifndef bmpg_uncc_edu_algorithms_Graph_Traits_hpp
#define bmpg_uncc_edu_algorithms_Graph_Traits_hpp

#include <string>
#include <limits>

namespace bmpg_uncc_edu {
namespace algorithms {
using namespace std;

/**
	This small file is a little tricky to document, so I (Hui Wang) am not
	going to bother. If you know it, you know it.
*/
template<bool C, typename T1, typename T2>
struct if_c
{
	typedef T1 type;	
	
};

template<typename T1, typename T2>
struct if_c<false,T1,T2>
{
	typedef T2 type;	
	
};



template<typename T>
struct is_primitive
{
	enum { value = false};
};


template<typename T>
struct is_primitive<T*>
{
	enum { value = true};
};


template<>
struct is_primitive<int>
{
	enum { value = true};
};


template<>
struct is_primitive<std::string>
{
	enum { value = true};
};


template<>
struct is_primitive<char>
{
	enum { value = true};
};



template<typename T,typename T2>
struct if_
{
private:
	typedef typename if_c<
			is_primitive<T>::value,
			T,
			T2
			>::type almost_type;
public:
	typedef almost_type type;
};

template<typename T>
struct Limit
{
public:
	T inf();
	T zero();
};


template<typename T>
T Limit<T>::inf()
{
	T t;
	return t.inf();
}

template<typename T>
T Limit<T>::zero()
{
	T t;
	return t.zero();
}

template<>
struct Limit<int>
{
public:
	int inf();
	int zero();
};



template<>
struct Limit<double>
{
public:
	double inf();
	double zero();
};


template<>
struct Limit<float>
{
public:
	float inf();
	float zero();
};

}	//namespace bmpg_uncc_edu::algorithms
}	//namespace bmpg_uncc_edu

#endif

