/**
   Copyright (c) 2008 by Hui Wang
   @author Hui Wang
   @date Sep 15, 2009
*/


#ifndef bmpg_uncc_edu_util_Enum2StringBinder_hpp
#define bmpg_uncc_edu_util_Enum2StringBinder_hpp

#include <sstream>
#include <boost/bimap.hpp>
#include <bmpg_uncc_edu/util/Exception.hpp>

namespace bmpg_uncc_edu {
namespace util {
using namespace std;

template<typename enum_t>
class Enum2StringBinder
{
public:
	typedef boost::bimap<enum_t, std::string> bimap_t;
	typedef typename bimap_t::value_type value_t;

	Enum2StringBinder();
	virtual void initialize() = 0;
	std::string get(const int& s) const;
	std::string get(const enum_t& s) const;
	enum_t get(const std::string& s) const;
	
	virtual ~Enum2StringBinder();
protected:
	void insert(const value_t& value);
private:
	bimap_t map_;
};


template<typename enum_t>
Enum2StringBinder<enum_t>::Enum2StringBinder()
{
}

template<typename enum_t>
enum_t Enum2StringBinder<enum_t>::get(const std::string& s) const
{
	typename bimap_t::right_map::const_iterator it;	
	it = map_.right.find(s);
	if(it == map_.right.end()){
		stringstream ss;
		ss << "Key \'" << s << "\' is not defined.";
		throw Exception(ss.str(), __FILE__, __LINE__);
	}
	return it->second;
}

template<typename enum_t>
std::string Enum2StringBinder<enum_t>::get(const int& s) const
{
	return get(static_cast<enum_t>(s));
}

template<typename enum_t>
std::string Enum2StringBinder<enum_t>::get(const enum_t& s) const
{
	typename bimap_t::left_map::const_iterator it;	
	it = map_.left.find(s);
	if(it == map_.left.end()){
		stringstream ss;
		ss << "Key " << s << " is not defined.";
		throw Exception(ss.str(), __FILE__, __LINE__);
	}
	return it->second;
}

template<typename enum_t>
void Enum2StringBinder<enum_t>::insert(const value_t& value)
{
	map_.insert(value);
}

template<typename enum_t>
Enum2StringBinder<enum_t>::~Enum2StringBinder()
{
	map_.clear();	
}


}	//namespace bmpg_uncc_edu::util
}	//namespace bmpg_uncc_edu

#endif

