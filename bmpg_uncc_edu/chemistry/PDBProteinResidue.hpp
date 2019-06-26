/**
   Copyright (c) 2008 by Mike Fairchild, Hui Wang
   @author: Mike Fairchild, http://www.mikef.org, mjfairch@gmail.com
   @author: Hui Wang
*/

#ifndef bmpg_uncc_edu_chemistry_PDBProteinResidue_hpp
#define bmpg_uncc_edu_chemistry_PDBProteinResidue_hpp

#include <vector>
#include <string>
#include <bmpg_uncc_edu/util/Handle.hpp>


namespace bmpg_uncc_edu {
	namespace chemistry {
		namespace library {
			class AminoAcid;
		}
	}
}

namespace bmpg_uncc_edu {
namespace chemistry {
using namespace std;
using namespace bmpg_uncc_edu::util;

class PDBProteinChain;	// forward declaration
class PDBProtein;	// forward declaration
class PDBAtom;

/**
 * PDBProteinResidue defines the details for a single residue in a protein, including all the atoms, energy, and structure. */


class PDBProteinResidue
{
public:
	typedef std::vector<PDBAtom *> atoms_t;					// Convenient typedef
	typedef atoms_t::iterator iterator;					// Convenient typedef
	typedef atoms_t::const_iterator const_iterator;				// Convenient typedef
	
	const bmpg_uncc_edu::chemistry::library::AminoAcid * aa;		// The amino acid type for this residue
	PDBProteinChain * chain;						// Parent PDBChain containing this residue
	PDBProtein * file;							// Parent PDBFile containing this residue
	Handle<PDBProteinResidue> next_res;					// The previous residue by number (NULL if non-existent)
	Handle<PDBProteinResidue> prev_res;					// The next residue by number (NULL if non-existent)
	Handle<PDBProteinResidue> sulf_res;					// The other residue forming the disulfide bond, if any.
	
	char chain_id;								// Together with number and insertion_code uniquely identifies a residue within a PDB file
	int number;								// Together with insertion_code uniquely identifies a residue within a chain
	char insertion_code;							// Together with number uniquely identifies a residue within a chain
	int user_id;								// User-defined id.  Define and use as you wish.
	
	double inter_cvdw_energy, inter_hbond_energy;
	
	PDBAtom * atom(const char * s,
		       bool try_alt_hydrogen_names = true) const;
	
	PDBAtom * atom(const char * symbol,
		       char remoteness,
		       char branch,
		       bool try_alt_hydrogen_names = true) const;
	
	PDBAtom * atom(const std::string & s,
		       bool try_alt_hydrogen_names = true) const
	{
		return atom(s.c_str(), try_alt_hydrogen_names);
	}
	
	PDBAtom * atom(const std::string & s,
		       char remoteness,
		       char branch,
		       bool try_alt_hydrogen_names = true) const
	{
		return atom(s.c_str(), remoteness, branch, try_alt_hydrogen_names);
	}
	
	bool has_atom(const std::string& s,
		      bool try_alt_hydrogen_names = true) const;
	
	size_t num_atoms() const { return atoms_.size(); }
	const atoms_t& atoms() const {return atoms_;}
	atoms_t& atoms() {return atoms_;}
	const_iterator begin() const{return atoms_.begin();}
	const_iterator end() const{return atoms_.end();}

	iterator begin() {return atoms_.begin();}
	iterator end() {return atoms_.end();}

	size_t hash_code();
	int nnc() const{return nnc_;}
	int& nnc(){return nnc_;}
	double asa() const {return asa_;}
	double& asa() {return asa_;}

	double& total_native_energy() {return total_native_energy_;}
	double total_native_energy() const {return total_native_energy_;}

	double& energy() {return total_native_energy_;}
	double energy() const {return total_native_energy_;}
	
	string& secondary_structure() {return sec_struct_;}
	string secondary_structure() const{return sec_struct_;}
	
	PDBProteinResidue();
	~PDBProteinResidue();

	
	friend std::ostream & operator<<(std::ostream & os, const PDBProteinResidue & res);
	bool operator<(const PDBProteinResidue & res) {return number < res.number;}
	
private:
	atoms_t atoms_;								// Set of PDBAtoms this residue contains
private:
	PDBProteinResidue(const PDBProteinResidue & rhs);			// hide copy constructor
	PDBProteinResidue & operator=(const PDBProteinResidue & rhs);		// hide assignment
	
	atoms_t atoms_with_all_occupancy(const char * s,
					 bool try_alt_hydrogen_names = false) const;
	
	int nnc_;								//number of neighboring contacts
	double asa_;								//accessible surface area, default to 0
	string sec_struct_;
	double total_native_energy_;
};


}	//namespace bmpg_uncc_edu::chemistry
}	//namespace bmpg_uncc_edu

#endif

