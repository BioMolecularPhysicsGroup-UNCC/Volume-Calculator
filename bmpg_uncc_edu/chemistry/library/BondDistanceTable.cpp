/**
   Copyright (c) 2008 by Hui Wang
   @author Hui Wang
   @date Jan 20, 2009
*/


#ifndef bmpg_uncc_edu_chemistry_library_BondDistanceTable_cpp
#define bmpg_uncc_edu_chemistry_library_BondDistanceTable_cpp

#include <cstdlib>
#include <iostream>
#include <vector>
#include <cmath>
#include <bmpg_uncc_edu/chemistry/Element.hpp>
#include <bmpg_uncc_edu/chemistry/library/BondDistanceTable.hpp>
#include <bmpg_uncc_edu/util/StringUtil.hpp>
#include <bmpg_uncc_edu/util/logger/LoggerFactory.hpp>
#include <bmpg_uncc_edu/chemistry/library/PeriodicTable.hpp>
#include <bmpg_uncc_edu/util/Exception.hpp>

namespace bmpg_uncc_edu {
namespace chemistry {
namespace library {
using namespace bmpg_uncc_edu::util;
using namespace bmpg_uncc_edu::util::logger;


BondDistanceTable* BondDistanceTable::instance_ = NULL;

BondDistanceTable::BondDistanceTable():loaded_(false)
{
}

BondDistanceTable* BondDistanceTable::instance()
{
	if(instance_ == NULL){
		instance_ = new BondDistanceTable();
	}
	return instance_;
}

/**
	Load the bond distance table from an input stream.
	@param input stream
*/
void BondDistanceTable::load(istream & is)
{
	if(loaded_)map_.clear();
	string line;
	char buf[1024];
	while(!is.eof()){
		is.getline(buf,1024);
		line = buf;
		if(line.empty())continue;
		string s1 = line.substr(0,2);
		string s2 = line.substr(2,2);
	
		PeriodicTable * table = PeriodicTable::instance();
		
		if(s1.at(1) == ' ')s1.erase(1);
		if(s2.at(1) == ' ')s2.erase(1);
		
		const Element* e1 = table->element(s1);
		const Element* e2 = table->element(s2);
		key_t key(e1,e2);
		
		string distance_string = line.substr(4);
	
		std::vector<string> dlist = StringUtil::tokenize(distance_string);
		std::vector<string>::iterator it;
		for(it = dlist.begin(); it != dlist.end(); it++){
			map_[key].push_back(atof((*it).c_str()));		//FIXME - what if multiple distances are provided?
		}
	}
	loaded_ = true;
}

/**
	Find the distance between two elements.
	@param Element 1
	@param Element 2
	@return iterator to the distance map
*/
double BondDistanceTable::find(const Element* e1,
			       const Element* e2,
			       size_t bond_number)
{
	if(!loaded_)
		throw bmpg_uncc_edu::util::Exception("BondDistanceTable::find(const Element* e1, const Element* e2): BondDistanceTable not loaded", __FILE__, __LINE__);
	
	
	BondDistanceTable::map_t::iterator it1, it2;
	it1 = map_.find(key_t(e1,e2));
	it2 = map_.find(key_t(e2,e1));
	if(it1 != map_.end()){			//switch and try again
		return (it1->second)[bond_number - 1];
	} else if(it2 != map_.end()){
		return (it2->second)[bond_number - 1];
	} else {
		throw bmpg_uncc_edu::util::Exception("BondDistanceTable::find(const Element* e1, const Element* e2): Elements not found in map", __FILE__, __LINE__);
	}
	
	//return it1;
}

/**
	Overloaded output stream operator <<.
*/
ostream& operator<<(ostream& out, const BondDistanceTable & rhs)
{
	BondDistanceTable::map_t::const_iterator it;
	for(it = rhs.map_.begin(); it != rhs.map_.end(); it++){
		const Element* e1 = it->first.first;
		const Element* e2 = it->first.second;
		if(e1 != NULL && e2 != NULL)
			out << e1->symbol() << " " << e2->symbol();
		BondDistanceTable::values_t::const_iterator itv;
		for(itv = (it->second).begin(); itv != (it->second).end(); itv++){
			out << " " << *itv;
		}
		out << "\n";
	}
	return out;
		
}

/**
	Free the memory allocated by the BondDistanceTable class.
*/
void BondDistanceTable::clear()
{
	map_.clear();
	loaded_ = false;
}
/**
	Destructor.
*/
BondDistanceTable::~BondDistanceTable()
{
	clear();
}


}	//namespace bmpg_uncc_edu::chemistry::library
}	//namespace bmpg_uncc_edu::chemistry
}	//namespace bmpg_uncc_edu

#endif

