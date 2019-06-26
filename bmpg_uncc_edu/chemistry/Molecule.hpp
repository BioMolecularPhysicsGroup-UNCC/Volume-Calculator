/**
   Copyright (c) 2008 by Hui Wang
   @author Hui Wang
   @date Jan 20, 2009
*/


#ifndef bmpg_uncc_edu_chemistry_Molecule_hpp
#define bmpg_uncc_edu_chemistry_Molecule_hpp

#include <vector>
#include <list>
#include <bmpg_uncc_edu/util/Handle.hpp>
#include <bmpg_uncc_edu/algorithms/Graph.hpp>

namespace bmpg_uncc_edu {
namespace chemistry {
using namespace bmpg_uncc_edu::util;
using namespace bmpg_uncc_edu::algorithms;

class Element;
class PDBAtom;
class CovalentBond;
class CovalentBonder;

/**
 * Molecule contains a list of PDBAtoms which are connected by covalent bonds. The molecule class owns the atoms. A 
 * user-defined covalent bonder is used to identify the covalent bonds among the atoms. The identified covalent bonds
 * are stored in the molecule.
 */
class Molecule
{
public:
	typedef std::vector<PDBAtom*> atoms_t;					//stores pointers
	typedef std::vector<const PDBAtom*> const_atoms_t;					//stores pointers
	typedef std::list<Handle<CovalentBond> > covalent_bonds_t;		//stores bond objects
	
	typedef Graph<const PDBAtom*,int> graph_t;
	typedef graph_t::vertex_list_t neighbor_list_t;
	
	int user_id;								// User defined id - define as you please
	
	Molecule();
	explicit Molecule(CovalentBonder* bonder);
	
	virtual void build_covalent_bonds();
	
	CovalentBonder* covalent_bonder() {return bonder_;}
	void covalent_bonder(CovalentBonder* bonder) {bonder_ = bonder;}
	const PDBAtom* closest_neighbor(const PDBAtom* atom,
					const Element* element,
					int degree) const;
	
	const PDBAtom* NthNearestNeighbor(const PDBAtom* atom_1,
					  const PDBAtom* atom_2,
					  int degree,
					  const PDBAtom* root_atom = NULL);
	
	const PDBAtom* NthNearestNeighbor(const PDBAtom* atom_1,
					  const Element* atom_2,
					  int degree,
					  const PDBAtom* root_atom = NULL);
	
	const PDBAtom* NthNearestNeighbor(const PDBAtom* atom_1,
					  const char* atom_2,
					  int degree,
					  const PDBAtom* root_atom = NULL);

	bool is_neighbor(const PDBAtom* atom1,
			 const PDBAtom* atom2) const;
	
	bool is_neighbor(const PDBAtom* atom,
			 const Element* element) const;
	
	size_t number_of_neighbors(const PDBAtom* atom) const;
	neighbor_list_t neighbors(const PDBAtom* atom) const;
	
	virtual atoms_t& atoms();
	
	covalent_bonds_t disulfide_bonds() const {return disulfide_bonds_;}
	
	size_t num_of_atoms() const;
	size_t num_of_covalent_bonds() const;
	bool has_atom(const PDBAtom* atom) const;
	bool has_bond(const PDBAtom* atom1,
		      const PDBAtom* atom2) const;
	
	double mass() const;							//total mass of the molecule
	
	const graph_t& graph() {return graph_;}
	
	void print_bonds(ostream & out);
	
	virtual ~Molecule();
protected:
	CovalentBonder* bonder_;						//covalent bonder
	atoms_t atoms_;								//vector of atoms, ordered by construction/addition
	covalent_bonds_t covalent_bonds_;					//list of covalent bonds
	covalent_bonds_t disulfide_bonds_;
	graph_t graph_;
};


}	//namespace bmpg_uncc_edu::chemistry
}	//namespace bmpg_uncc_edu

#endif

