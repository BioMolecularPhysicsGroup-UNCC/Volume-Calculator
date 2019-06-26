/**
   Copyright (c) 2008 by Hui Wang
   @author Hui Wang
   @date Jan 20, 2009
*/


#ifndef bmpg_uncc_edu_chemistry_Molecule_cpp
#define bmpg_uncc_edu_chemistry_Molecule_cpp

#include <bmpg_uncc_edu/chemistry/PDBAtom.hpp>
#include <bmpg_uncc_edu/chemistry/CovalentBond.hpp>
#include <bmpg_uncc_edu/chemistry/CovalentBonder.hpp>
#include <bmpg_uncc_edu/chemistry/Molecule.hpp>
#include <bmpg_uncc_edu/chemistry/helper/PDBHelper.hpp>
#include <bmpg_uncc_edu/util/Exception.hpp>

namespace bmpg_uncc_edu {
namespace chemistry {
using namespace bmpg_uncc_edu::chemistry::helper;

/**
	Default constructor. A covalent bonder whoes hash grid is initialized with grid size 5 angstrom is created.
*/
Molecule::Molecule() :
		   user_id(0)
{
	bonder_ = new CovalentBonder(5);
}

/**
	Constructor.
	@param bonder CovalentBonder
*/
Molecule::Molecule(CovalentBonder* bonder) :
		   bonder_(bonder)
{

}
	
	
/**
	Destructor. A molecule frees the memory of the atoms.
*/
Molecule::~Molecule()
{	
	atoms_t::iterator aiter;
	for (aiter = atoms_.begin(); aiter != atoms_.end(); ++aiter){
		delete *aiter;
	}
	atoms_.clear();
	
	covalent_bonds_.clear();
	
	if(bonder_ != NULL)
		delete bonder_, bonder_ = NULL;
}

/**
	Build covalent bonds using the covalent bonder and store in the class.
*/
void Molecule::build_covalent_bonds()
{
	if(bonder_ == NULL)
		throw bmpg_uncc_edu::util::Exception("Molecule::build_covalent_bonds(): set covalent_bonder first.", __FILE__, __LINE__);
	
	covalent_bonds_.clear();
	
	covalent_bonds_ = bonder_->build_covalent_bonds(*this);
	disulfide_bonds_ = bonder_->disulfide_bonds();
	
	atoms_t::iterator ita;
	for(ita = atoms_.begin(); ita != atoms_.end(); ita++){
		graph_.add_vertex(*ita);	
	}
	
	covalent_bonds_t::iterator it;
	for(it = covalent_bonds_.begin(); it != covalent_bonds_.end(); it++){
		Handle<CovalentBond> bond = *it;	
		graph_.add_edge(bond->first_atom(), bond->second_atom());
	}
}

const PDBAtom* Molecule::NthNearestNeighbor(const PDBAtom* atom_1,
					    const PDBAtom* atom_2,
					    int degree,
					    const PDBAtom* root_atom)
{
	const PDBAtom* atom = NULL;

	if( degree == 0 )
		return NULL;
	neighbor_list_t neighbor_list = graph_.neighbors(atom_1);
	
	neighbor_list_t::iterator next_atom = neighbor_list.begin();
	while(next_atom != neighbor_list.end()){
		if(*next_atom == atom_2)
			return(*next_atom);
		else if(*next_atom == root_atom)
			next_atom++;
		else{
			atom = NthNearestNeighbor(*next_atom, atom_2, degree-1, atom_1);
			if(atom != NULL)
				return(atom);
			else
				next_atom++;
		}
	}

	return NULL;
}


const PDBAtom* Molecule::NthNearestNeighbor(const PDBAtom* atom_1,
					    const Element* atom_2,
					    int degree,
					    const PDBAtom* root_atom)
{
	const PDBAtom* atom = NULL;

	if( degree == 0 )
		return NULL;
	
	neighbor_list_t neighbor_list = graph_.neighbors(atom_1);
	neighbor_list_t::iterator next_atom = neighbor_list.begin();
	while(next_atom != neighbor_list.end()){
		if(PDBHelper::atom_to_element(*next_atom) == atom_2 )
			return(*next_atom);
		else if(*next_atom == root_atom)
			next_atom++;
		else{
			atom = NthNearestNeighbor(*next_atom, atom_2, degree-1, atom_1);
			if(atom != NULL)
				return(atom);
			else
				next_atom++;
		}
	}

	return NULL;
}


const PDBAtom* Molecule::NthNearestNeighbor(const PDBAtom* atom_1,
					    const char* atom_2,
					    int degree,
					    const PDBAtom* root_atom)
{
	const PDBAtom* atom = NULL;

	if( degree == 0 )
		return NULL;
	string atom_name_2 = atom_2;
	
	neighbor_list_t neighbor_list = graph_.neighbors(atom_1);
	
	neighbor_list_t::iterator next_atom = neighbor_list.begin();
	while(next_atom != neighbor_list.end()){
		string name1 = (*next_atom)->atom_name;
		if(name1 == atom_name_2)
			return(*next_atom);
		else if(*next_atom == root_atom)
			next_atom++;
		else{
			atom = NthNearestNeighbor(*next_atom, atom_2, degree-1, atom_1);
			if(atom != NULL)
				return(atom);
			else
				next_atom++;
		}
	}

	return NULL;
}


/**
	Test if two atoms are neighbors based on the record in the library.
	@param atom1 Atom 1
	@param atom2 Atom 2
	@return true if neighbors; false if not
*/
bool Molecule::is_neighbor(const PDBAtom* atom1,
			   const PDBAtom* atom2) const
{
	return graph_.is_neighbor(atom1,atom2);
}


Molecule::neighbor_list_t Molecule::neighbors(const PDBAtom* atom) const
{
	return graph_.neighbors(atom);
}

/**
	Test if an atom is neighboring with an element.
	@param atom atom anme
	@param element element
	@return true if neighbors; false if not
*/
bool Molecule::is_neighbor(const PDBAtom* atom,
			   const Element* element) const
{
	neighbor_list_t nb_list = graph_.neighbors(atom);
	
	neighbor_list_t::const_iterator it;
	for(it = nb_list.begin(); it != nb_list.end(); it++){
		const PDBAtom* nbr = *it;
		if(nbr->element() == element)
			return true;
	}
	return false;
}

size_t Molecule::number_of_neighbors(const PDBAtom* atom) const
{
	return graph_.number_of_neighbors(atom);	
}

/**
	Returns list of all the atoms in the molecule.
	@return list of atoms
*/
Molecule::atoms_t& Molecule::atoms()
{
	return atoms_;
}

/**
	Returns the total number of atoms in the molecule.
	@return total number of atoms
*/
size_t Molecule::num_of_atoms() const
{
	return atoms_.size();
}

/**
	Returns the total number of covalent bonds in the molecule.
	@return total number of covalent bonds
*/
size_t Molecule::num_of_covalent_bonds() const
{
	return covalent_bonds_.size();
}

/**
	Check if the molecule has a given atom.
	@param aomt Atom pointer
	@return true if atom is found; false if not
*/
bool Molecule::has_atom(const PDBAtom* atom) const
{
	atoms_t::const_iterator iter;
	for(iter = atoms_.begin(); iter != atoms_.end(); iter++){
		if(atom == *iter)return true;
	}
	return false;
}

/**
	Check if the molecule contains a covalent bond between two given atoms.
	@param atom1 Atom 1
	@param atom2 Atom 2
	@return true if there is a bond between the two atoms; false if not
*/
bool Molecule::has_bond(const PDBAtom* atom1,
			const PDBAtom* atom2) const
{
	covalent_bonds_t::const_iterator it;
	for(it = covalent_bonds_.begin(); it != covalent_bonds_.end(); it++){
		const PDBAtom* a1 = (*it)->first_atom();
		const PDBAtom* a2 = (*it)->second_atom();
		if( (a1 == atom1 && a2 == atom2) || (a1 == atom2 && a2 == atom1))
                    return true;
	}
	return false;
}

/**
	Returns the total mass of the molecule (in units of atomic mass u).
	@return total mass in units of atomic mass u
*/
double Molecule::mass() const
{	
	atoms_t::const_iterator iter;
	double m = 0;
	for(iter = atoms_.begin(); iter != atoms_.end(); iter++){
		m += (*iter)->mass();
	}
	return m;
}

/**
	Print the covalent bonds in the molecule to an output stream.
	@param out output stream
*/
void Molecule::print_bonds(ostream & out)
{
	covalent_bonds_t::iterator it;
	for(it = covalent_bonds_.begin(); it != covalent_bonds_.end(); it++){
		out << **it << endl;
	}
}


}	//namespace bmpg_uncc_edu::chemistry
}	//namespace bmpg_uncc_edu

#endif

