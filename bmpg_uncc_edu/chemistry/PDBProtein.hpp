/**
   Copyright (c) 2008 by Mike Fairchild, Hui Wang
   @author Mike Fairchild, http://www.mikef.org, mjfairch@gmail.com
   @author Hui Wang
*/

#ifndef bmpg_uncc_edu_chemistry_PDBProtein_hpp
#define bmpg_uncc_edu_chemistry_PDBProtein_hpp

#include <string>
#include <bmpg_uncc_edu/util/Handle.hpp>
#include <bmpg_uncc_edu/chemistry/PDBDefs.hpp>
#include <bmpg_uncc_edu/chemistry/PDBProteinChain.hpp>
#include <bmpg_uncc_edu/chemistry/Complex.hpp>
#include <bmpg_uncc_edu/chemistry/detail/PDBChainManager.hpp>


namespace bmpg_uncc_edu {
	namespace fast {
		namespace input {
			class EnvironmentCondition;
		}
	}
}

namespace bmpg_uncc_edu {
namespace chemistry {
using namespace std;
using namespace bmpg_uncc_edu::fast::input;

class PDBAtom;
class PDBProteinResidue;

/**
 * PDBProtein reads a protein from a PDB file, and stores and calculates all the detailed information found in this file.
 * The PDBProtein is constructed one Atom class at a time, also building the covalent bonds and tagging the beginning and end
 * of each residue.  In this way, all the topological information is accurately captured for the protein.
 */

class PDBProtein : public Complex<PDBProteinChain>
{
public:
	typedef Molecule::covalent_bonds_t covalent_bonds_t;
	
	static const char all_chains;
	
	std::string filename;					// User-defined.  Should contain the PDB filename.
	int user_id;						// User-defined.  Define as you please.
	PDBProtein();
	~PDBProtein();

	size_t num_atoms() const;

	Handle<EnvironmentCondition> environment_condition() const {return condition_;}
//	Handle<EnvironmentCondition>& set_environment_condition(Handle<EnvironmentCondition> cond);
	
	class ResidueIterator;
	ResidueIterator residue_begin() const;
	ResidueIterator residue_end() const;
	ResidueIterator residue_begin();
	ResidueIterator residue_end();
	typedef ResidueIterator residue_iterator_t;
	
	const PDBAtom * atom(int n) const;
	const PDBAtom * atom_by_user_id(int id) const;
	Handle<PDBProteinResidue> residue(int n, char chain_id = PDB_DEFAULT_CHAIN_ID, char insertion_code = PDB_DEFAULT_INSERTION_CODE) const;
	Handle<PDBProteinResidue> residue_by_user_id(int id) const;
	
	void build_covalent_bonds();
	bool has_default_chain() const;
	bool has_only_default_chain() const;
	PDBProteinChain * default_chain() const;
	void keep_chain_only(const char& id);

	
	void remove_empty_chains();
	molecules_map_t::iterator get_empty_chain();
	
	covalent_bonds_t disulfide_bonds() const;

	size_t num_molecules() const;
	size_t num_residues() const;
	bool discontiguous_residue_found() const {return discontiguous_residue_found_;}
	
	void read(const string& fname);
	void read(std::istream & is);
	void write(std::ostream & os) const;
	friend std::ostream & operator<<(std::ostream & os, const PDBProtein & file);
	friend std::istream & operator>>(std::istream & is, PDBProtein & file);

protected:
	PDBAtom* read_atom(const std::string & line);
	void write_atom(const std::ostream & os);
	void crosslink(PDBAtom* atm, Handle<PDBProteinResidue> &pres);

private:
	PDBProtein(const PDBProtein & rhs); 					// Hide copy constructor
	PDBProtein & operator=(const PDBProtein & rhs);  			// Hide assignment operator
	
	bool discontiguous_residue_found_, decrease_residue_number_found;
	Handle<detail::PDBChainManager> chain_manager_;
	char current_chain_;
	Handle<EnvironmentCondition> condition_;
};

class PDBProtein::ResidueIterator
{
public:
	ResidueIterator();
	ResidueIterator(const PDBProtein* protein);
	ResidueIterator(const PDBProtein* protein, molecules_map_t::const_iterator itc);
	ResidueIterator(const PDBProtein* protein, molecules_map_t::const_iterator itc, residues_t::const_iterator itr);
	ResidueIterator& operator++();						//postfix operator++
	ResidueIterator& operator++(int);					//a dummy int overloads prefix operator++	
	bool operator==(const PDBProtein::ResidueIterator& rhs);	
	bool operator!=(const PDBProtein::ResidueIterator& rhs);	
	ResidueIterator& operator=(const PDBProtein::ResidueIterator& rhs);	
	Handle<PDBProteinResidue> operator*();	
private:
	const PDBProtein* protein_;
	molecules_map_t::const_iterator curr_molecule;
	residues_t::const_iterator curr_res;
	
};


}	//namespace bmpg_uncc_edu::chemistry
}	//namespace bmpg_uncc_edu

#endif

