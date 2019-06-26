/**
   Copyright (c) 2008 by Hui Wang
   @author Hui Wang
   @date Jan 20, 2009
*/


#ifndef bmpg_uncc_edu_chemistry_CovalentBond_hpp
#define bmpg_uncc_edu_chemistry_CovalentBond_hpp

#include <bmpg_uncc_edu/chemistry/Bond.hpp>

namespace bmpg_uncc_edu {
namespace chemistry {

class PDBAtom;

/**
 * CovalentBond contains Bond information between two PDBAtoms specific to a covalent bond. 
 */

class CovalentBond : public Bond
{
public:
	CovalentBond(const PDBAtom* atom1,
		     	 const PDBAtom* atom2);
	
	const PDBAtom* first_atom() const;
	const PDBAtom* second_atom() const;
	
	Bond::bond_t bond_type() {return Bond::COVALENT;}
	size_t num_of_atoms() const {return 2;}
	
	friend ostream& operator<<(ostream & out, const CovalentBond & rhs);
private:
	const PDBAtom *first_atom_, *second_atom_;
};


}	//namespace bmpg_uncc_edu::chemistry
}	//namespace bmpg_uncc_edu

#endif

