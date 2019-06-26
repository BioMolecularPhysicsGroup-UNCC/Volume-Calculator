/**
   Copyright (c) 2008 by Hui Wang and Brandon Hespenheide 
   @author Hui Wang
   @date Jan 20, 2009
*/


#ifndef bmpg_uncc_edu_chemistry_CovalentBonder_hpp
#define bmpg_uncc_edu_chemistry_CovalentBonder_hpp

#include <limits>
#include <list>
#include <bmpg_uncc_edu/util/Handle.hpp>

namespace bmpg_uncc_edu {
namespace chemistry {
using namespace bmpg_uncc_edu::util;

class PDBAtom;
class PDBProtein;
class Molecule;
class CovalentBond;

/**
 * CovalentBonder identifies the covalent bonds of a molecule. Once the bonds are created, CovalentBonder does NOT store the bonds.
 * The calling object is responsible for deallocating the memory of the bonds.
*/
class CovalentBonder
{
public:
	typedef std::list<Handle<CovalentBond> > cov_bonds_t;
	typedef std::list<Handle<CovalentBond> > disulfide_bonds_t;

	explicit CovalentBonder(double l) :
		       				grid_length_(l),
		       				built_(false)
		       				{
		       				}
	
	virtual cov_bonds_t build_covalent_bonds(Molecule & mol);
	virtual cov_bonds_t build_covalent_bonds(PDBProtein& protein);
	disulfide_bonds_t disulfide_bonds() const;
	
	bool meets_atom_atom_distance_cutoff(const PDBAtom* atom1,
					     const PDBAtom* atom2);
	
	double& hard_distance_cutoff() {return hard_distance_cutoff_;}
	double hard_distance_cutoff() const {return hard_distance_cutoff_;}
	
	bool check_peptide_bond(const PDBAtom * atom1,
				const PDBAtom * atom2);
	
	bool check_disulfide_bond(const PDBAtom * atom1,
				  const PDBAtom * atom2);
	
	bool check_internucleic_bond(const PDBAtom * atom1,
				     const PDBAtom * atom2);
	
	virtual ~CovalentBonder(){}
private:
	double grid_length_;
	static double hard_distance_cutoff_;
	disulfide_bonds_t disulfide_bonds_;
	bool built_;
	
};

}	//namespace bmpg_uncc_edu::chemistry
}	//namespace bmpg_uncc_edu

#endif

