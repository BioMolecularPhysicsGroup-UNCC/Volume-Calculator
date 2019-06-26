/**
   Copyright (c) 2008 by Hui Wang
   @author Hui Wang
   @date Jan 20, 2009
*/


#ifndef bmpg_uncc_edu_chemistry_Complex_hpp
#define bmpg_uncc_edu_chemistry_Complex_hpp

#include <list>
#include <map>
#include <vector>
#include <bmpg_uncc_edu/util/Handle.hpp>

namespace bmpg_uncc_edu {
namespace chemistry {
using namespace std;
using namespace bmpg_uncc_edu::util;

namespace hbond {
	class HydrogenBond;
}

class PDBAtom;
class PDBProteinResidue;
class ResiduePair;
class Ion;

template <class Unit>

/**
 * Complex contains one or more covalently bonded Molecules.
*/

class Complex
{
public:
	typedef Unit molecule_t;
	typedef typename std::vector<Unit*> molecules_t;
	typedef std::vector<PDBAtom*> atoms_t;
	typedef std::vector<Handle<PDBProteinResidue> > residues_t;		// Convenient typedef
	
	typedef std::vector<Ion*> ions_t;
	typedef typename std::map<char, Unit*> molecules_map_t;			//Convenient typedef
	typedef std::vector<Handle<hbond::HydrogenBond> > hydrogen_bonds_t;		//stores bond objects
	typedef std::list<Handle<ResiduePair> > residue_pairs_t;
	
	class AtomIterator;
	AtomIterator atom_begin() const;
	AtomIterator atom_end() const;
	AtomIterator atom_begin();
	AtomIterator atom_end();

	typedef AtomIterator atom_iterator_t;
	
	typedef typename molecules_map_t::iterator molecule_iterator_t;
	typedef typename molecules_map_t::const_iterator molecule_const_iterator_t;
	
	molecule_iterator_t molecule_begin() {return molecules_.begin();}
	molecule_iterator_t molecule_end() {return molecules_.end();}
	molecule_const_iterator_t molecule_begin() const{return molecules_.begin();}
	molecule_const_iterator_t molecule_end() const{return molecules_.end();}
	
	
	
	molecules_t molecules() const;
	Unit* molecule(char chain_id) const;
	Unit* molecule_by_user_id(int id) const;
	
	void assign_hydrogen_bonds(const hydrogen_bonds_t& hbonds);
	hydrogen_bonds_t hydrogen_bonds() const {return hydrogen_bonds_;}
	
	void assign_residue_pairs(const residue_pairs_t& res_pairs);
	residue_pairs_t residue_pairs() const {return residue_pairs_;}
	
	int id() const {return id_;}
	
	atoms_t atoms() const;
	size_t num_atoms() const;
	size_t num_molecules() const { return molecules_.size(); }
	
	Complex();
	~Complex();
	
protected:
	molecules_map_t molecules_;						//Set of PDBProteinChains for this PDB file
	
private:
	int id_;
	atoms_t atoms_;
	ions_t ions_;
	hydrogen_bonds_t hydrogen_bonds_;
	residue_pairs_t residue_pairs_;
};


template <class Unit>
class Complex<Unit>::AtomIterator
{
public:
	AtomIterator();
	AtomIterator(const Complex<Unit>* complex);
	AtomIterator(const Complex<Unit>* complex,
		     typename molecules_map_t::const_iterator itc);
	
	AtomIterator(const Complex<Unit>* complex,
		     typename molecules_map_t::const_iterator itc,
		     typename atoms_t::const_iterator itr);
	
	AtomIterator& operator++();						//postfix operator++
	AtomIterator& operator++(int);						//a dummy int overloads prefix operator++
	
	bool operator==(const typename Complex<Unit>::AtomIterator& rhs);	
	bool operator!=(const typename Complex<Unit>::AtomIterator& rhs);	
	AtomIterator& operator=(const typename Complex<Unit>::AtomIterator& rhs);	
	
	PDBAtom* operator*();
	
private:
	const Complex<Unit>* complex_;
	typename molecules_map_t::const_iterator curr_molecule;
	typename atoms_t::const_iterator curr_atom;
};

}	//namespace bmpg_uncc_edu::chemistry
}	//namespace bmpg_uncc_edu

#endif

