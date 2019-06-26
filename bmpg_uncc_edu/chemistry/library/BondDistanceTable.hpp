/**
   Copyright (c) 2008 by Hui Wang
   @author Hui Wang
   @date Jan 20, 2009
*/


#ifndef bmpg_uncc_edu_chemistry_library_BondDistanceTable_hpp
#define bmpg_uncc_edu_chemistry_library_BondDistanceTable_hpp

#include <map>
#include <vector>

namespace bmpg_uncc_edu {
	namespace chemistry {
		class Element;
	}
}

namespace bmpg_uncc_edu {
namespace chemistry {
namespace library {
using namespace std;

/**
 * BondDistanceTable defines the standard distances between atoms that are covalently bonded.
 * PeriodicTable must be loaded before loading BondDistanceTable.
*/
class BondDistanceTable
{
public:
	typedef std::pair<const Element*, const Element*> key_t;
	typedef std::vector<double> values_t;
	typedef std::map<key_t, values_t> map_t;

	static BondDistanceTable* instance();
	
	void load(istream & is);
	double find(const Element* e1, 
		    const Element* e2,
		    size_t bond_number = 1);	//number of bonds connecting the two elements
	
	map_t::iterator end() {return map_.end();}
	friend ostream& operator<<(ostream& out, const BondDistanceTable & rhs);
	void clear();
	
	virtual ~BondDistanceTable();	
private:
	BondDistanceTable();
	BondDistanceTable(const BondDistanceTable & rhs);
	BondDistanceTable& operator=(const BondDistanceTable & rhs);
private:
	map_t map_;
	static BondDistanceTable * instance_;
	bool loaded_;
};



}	//namespace bmpg_uncc_edu::chemistry::library
}	//namespace bmpg_uncc_edu::chemistry
}	//namespace bmpg_uncc_edu

#endif

