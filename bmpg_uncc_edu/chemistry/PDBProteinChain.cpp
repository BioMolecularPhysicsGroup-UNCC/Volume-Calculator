/**
   Copyright (c) 2008 by Mike Fairchild, Hui Wang
   @author: Mike Fairchild, http://www.mikef.org, mjfairch@gmail.com
   @author: Hui Wang
*/

#ifndef bmpg_uncc_edu_chemistry_PDBProteinChain_cpp
#define bmpg_uncc_edu_chemistry_PDBProteinChain_cpp

#include <bmpg_uncc_edu/chemistry/PDBDefs.hpp>
#include <bmpg_uncc_edu/chemistry/PDBAtom.hpp>
#include <bmpg_uncc_edu/chemistry/PDBProtein.hpp>
#include <bmpg_uncc_edu/chemistry/PDBProteinResidue.hpp>
#include <bmpg_uncc_edu/chemistry/PDBProteinChain.hpp>

namespace bmpg_uncc_edu {
namespace chemistry {
using namespace std;

PDBProteinChain::PDBProteinChain() : file(NULL), chain_id(PDB_DEFAULT_CHAIN_ID)
{
	// Reserve space for an "average" number of residues in a chain.
	// I admit this is a guess.  The datastructure will grow
	// automatically as needed.  This call to reserve() improves performance
	// at the beginning.
	residues.reserve(50);
}

PDBProteinChain::~PDBProteinChain()
{
	residues.clear();
}

//implement atoms() function in Molecule
PDBProteinChain::atoms_t& PDBProteinChain::atoms()
{
	return atoms_;
}


void PDBProteinChain::remove_atom(const PDBAtom* atm)
{
	if(atm == NULL || !has_atom(atm))
		return;

	atoms_t::iterator it, it_found;
	for(it = atoms_.begin(); it != atoms_.end(); it++){			//FIXME - XXXX
		if(atm == *it){
			it_found = it;
			break;
		}
	}
	
	
	if(it_found != atoms_.end()){
		atoms_.erase(it_found);
	}
	delete atm, atm = NULL;
}


bool PDBProteinChain::remove_residue(Handle<PDBProteinResidue> res)
{
	if(res.empty() || !has_residue(res))
		return false;
	
	Handle<PDBProteinResidue> null_res;
	
	if(!(res->prev_res).empty()){
		res->prev_res->next_res = null_res;
	}
	
	
	if(!(res->next_res).empty()){
		res->next_res->prev_res = null_res;
	
	}
	
	atoms_t::iterator iter;
	for(iter = res->atoms().begin(); iter != res->atoms().end(); iter++){
		remove_atom(*iter);
	}
	
	residues_t::iterator it, it_found;
	for(it = residues.begin(); it != residues.end(); it++){
		if(res == *it){
			it_found = it;
			break;
		}
	}
	if(it_found != residues.end())
		residues.erase(it_found);
        
	return true;
}


bool PDBProteinChain::has_atom(const PDBAtom* atm) const
{
	bool found = false;
	atoms_t::const_iterator iter;
	for (iter = atoms_.begin(); iter != atoms_.end(); ++iter){
		if(atm == *iter)
			found = true;
	}
	return found;
}

bool PDBProteinChain::has_residue(Handle<PDBProteinResidue> res) const
{
	bool found = false;
	residues_t::const_iterator riter;
	for (riter = residues.begin(); riter != residues.end(); ++riter){
		if(res == *riter)
			found = true;
	}
	return found;
}


}	//namespace bmpg_uncc_edu::chemistry
}	//namespace bmpg_uncc_edu

#endif

