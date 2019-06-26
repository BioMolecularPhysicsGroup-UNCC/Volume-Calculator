/**
   Copyright (c) 2008 by Hui Wang
   @author Hui Wang
   @date Jan 20, 2009
*/


#ifndef bmpg_uncc_edu_chemistry_Bond_hpp
#define bmpg_uncc_edu_chemistry_Bond_hpp

#include <iosfwd>

namespace bmpg_uncc_edu {
namespace chemistry {
using namespace std;

/**
 * Bond is the base class for either a hydrogen or covalent bond
 */

class Bond
{
public:
	typedef enum {HYDROGEN, COVALENT} bond_t;
	virtual bond_t bond_type() = 0;
	virtual size_t num_of_atoms() const = 0;
	virtual ~Bond(){}
protected:	
};

}	//namespace bmpg_uncc_edu::chemistry
}	//namespace bmpg_uncc_edu

#endif

