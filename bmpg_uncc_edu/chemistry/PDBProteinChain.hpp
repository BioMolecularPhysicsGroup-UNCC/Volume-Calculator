/**
   Copyright (c) 2008 by Mike Fairchild, Hui Wang
   @author: Mike Fairchild, http://www.mikef.org, mjfairch@gmail.com
   @author: Hui Wang
*/

#ifndef bmpg_uncc_edu_chemistry_PDBProteinChain_hpp
#define bmpg_uncc_edu_chemistry_PDBProteinChain_hpp

#include <vector>
#include <bmpg_uncc_edu/util/Handle.hpp>
#include <bmpg_uncc_edu/chemistry/PDBDefs.hpp>
#include <bmpg_uncc_edu/chemistry/Molecule.hpp>


namespace bmpg_uncc_edu {
namespace chemistry {
using namespace std;

class PDBProtein; 								// forward declaration
class PDBAtom;
class PDBProteinResidue;

/**
 * PDBProteinChain is a protein-specific implementation of Molecule, defining the covalently bonded residues of a protein.
 */

class PDBProteinChain : public Molecule
{
public:
	typedef Handle<PDBProteinResidue> residue_t;                            // Convenient typedef
	typedef std::vector<Handle<PDBProteinResidue> > residues_t;		// Convenient typedef
	typedef std::vector<PDBAtom *> atoms_t;					// Convenient typedef
	typedef residues_t::const_iterator residue_const_iterator_t;
	typedef residues_t::const_iterator residue_iterator_t;
	
	
	PDBProtein * file;							// Parent PDBFile to which this chain belongs
	char chain_id;							
	residues_t residues;							// Collection of residues for this chain
	
	bool is_default() const { return (chain_id == PDB_DEFAULT_CHAIN_ID); }
	atoms_t& atoms();							//implement atoms() function in Molecule
	
	residue_iterator_t residue_begin() {return residues.begin();}
	residue_iterator_t residue_end() {return residues.end();}
	residue_const_iterator_t residue_begin() const{return residues.begin();}
	residue_const_iterator_t residue_end() const{return residues.end();}
	
	bool remove_residue(Handle<PDBProteinResidue> res);

	bool has_residue(Handle<PDBProteinResidue> res) const;
	bool has_atom(const PDBAtom* atm) const;

	PDBProteinChain();
	~PDBProteinChain();
	
private:
	void remove_atom(const PDBAtom* atm);						//can only be called by remove_residue()
	
	PDBProteinChain(const PDBProteinChain & rhs);				// hide copy constructor
	PDBProteinChain & operator=(const PDBProteinChain & rhs);		// hide assignment
};
	
}	//namespace bmpg_uncc_edu::chemistry
}	//namespace bmpg_uncc_edu

#endif

