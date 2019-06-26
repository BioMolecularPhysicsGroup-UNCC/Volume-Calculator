/**
   Copyright (c) 2008 by Hui Wang
   @author Hui Wang
   @date Jan 20, 2009
*/


#ifndef bmpg_uncc_edu_chemistry_Complex_cpp
#define bmpg_uncc_edu_chemistry_Complex_cpp

#include <iostream>
#include <bmpg_uncc_edu/chemistry/Complex.hpp>
#include <bmpg_uncc_edu/chemistry/Molecule.hpp>
#include <bmpg_uncc_edu/chemistry/PDBProteinChain.hpp>
#include <bmpg_uncc_edu/algorithms/HashGrid.hpp>
#include <bmpg_uncc_edu/util/SequentialIdGenerator.hpp>
#include <bmpg_uncc_edu/util/logger/LoggerFactory.hpp>

namespace bmpg_uncc_edu {
namespace chemistry {
using namespace std;
using namespace bmpg_uncc_edu::algorithms;
using namespace bmpg_uncc_edu::util;
using namespace bmpg_uncc_edu::util::logger;

template <class Unit>
Complex<Unit>::Complex() :
		       id_(SequentialIdGenerator::get())
{
	
}
	
template <class Unit>
Complex<Unit>::~Complex()
{
	// Delete memory we allocated and clear out data structures
	typename molecules_map_t::iterator citer;
	for (citer = molecules_.begin(); citer != molecules_.end(); ++citer){
		delete citer->second;
	}
	molecules_.clear();
}

template <class Unit>
typename Complex<Unit>::molecules_t Complex<Unit>::molecules() const
{
	molecules_t result;
	typename molecules_map_t::const_iterator iter;
	for (iter = molecules_.begin(); iter != molecules_.end(); ++iter) {
		result.push_back(iter->second);
	}
	return result;
}

template <class Unit>
Unit* Complex<Unit>::molecule(char chain_id) const
{
	typename molecules_map_t::const_iterator iter;
	for (iter = molecules_.begin(); iter != molecules_.end(); ++iter) {
		if (iter->first == chain_id){
			return iter->second;
		}
	}
	
	return NULL;
}

template <class Unit>
void Complex<Unit>::assign_hydrogen_bonds(const hydrogen_bonds_t& hbonds)
{
	hydrogen_bonds_ = hbonds;
}

template <class Unit>
void Complex<Unit>::assign_residue_pairs(const residue_pairs_t& res_pairs)
{
	residue_pairs_ = res_pairs;
}
	

template <class Unit>
Unit* Complex<Unit>::molecule_by_user_id(int id) const
{
	typename molecules_map_t::const_iterator iter;
	for (iter = molecules_.begin(); iter != molecules_.end(); ++iter) {
		if (iter->second->user_id == id)
			return iter->second;
	}
	return NULL;
}


template <class Unit>
typename Complex<Unit>::atoms_t Complex<Unit>::atoms() const
{
	atoms_t result;
	atom_iterator_t iter;
	for(iter = atom_begin(); iter != atom_end(); iter++){
		result.push_back(*iter);
	}
	return result;
}


template <class Unit>
size_t Complex<Unit>::num_atoms() const
{
	size_t total = 0;
	typename molecules_map_t::const_iterator it;
	for(it = molecules_.begin(); it != molecules_.end(); it++){
		total += it->second->atoms().size();
	}
	return total;
}

template <class Unit>
Complex<Unit>::AtomIterator::AtomIterator() : complex_(NULL)
{
}

template <class Unit>
Complex<Unit>::AtomIterator::AtomIterator(const Complex* p) :
					  complex_(p)
{
	if(complex_->molecules_.size() == 0){
		curr_molecule = complex_->molecules_.end();
		return;
	}
	curr_molecule = complex_->molecules_.begin();
	
	while (curr_molecule != complex_->molecules_.end() && (curr_molecule->second == NULL || curr_molecule->second->atoms().size() == 0)){
			curr_molecule++;
	}

	if(curr_molecule != complex_->molecules_.end())
		curr_atom = curr_molecule->second->atoms().begin();
}

template <class Unit>
Complex<Unit>::AtomIterator::AtomIterator(const Complex<Unit>* p,
					  typename molecules_map_t::const_iterator it) :
					  complex_(p),
					  curr_molecule(it)
{
}

template <class Unit>
Complex<Unit>::AtomIterator::AtomIterator(const Complex<Unit>* p,
					  typename molecules_map_t::const_iterator it,
					  typename atoms_t::const_iterator itr) :
					  complex_(p),
					  curr_molecule(it),
					  curr_atom(itr)
{
}

template <class Unit>
typename Complex<Unit>::AtomIterator Complex<Unit>::atom_begin() const
{
	return AtomIterator(this);
}

template <class Unit>
typename Complex<Unit>::AtomIterator Complex<Unit>::atom_end() const
{
	return AtomIterator(this,molecules_.end());
}

template <class Unit>
typename Complex<Unit>::AtomIterator Complex<Unit>::atom_begin()
{
	return AtomIterator(this);
}


template <class Unit>
typename Complex<Unit>::AtomIterator Complex<Unit>::atom_end()
{
	return AtomIterator(this,molecules_.end());
}

template <class Unit>
PDBAtom* Complex<Unit>::AtomIterator::operator*()
{
	return *curr_atom;
}

template <class Unit>
typename Complex<Unit>::AtomIterator& Complex<Unit>::AtomIterator::operator++()
{		
	if(curr_atom != curr_molecule->second->atoms().end())
		curr_atom++;

	if(curr_atom == curr_molecule->second->atoms().end()){
		curr_molecule++;
		while (curr_molecule != complex_->molecules_.end() && curr_molecule->second->atoms().size() == 0){
			curr_molecule++;
		}

		if(curr_molecule != complex_->molecules_.end())
			curr_atom = curr_molecule->second->atoms().begin();
	}

	return *this;
}

template <class Unit>
typename Complex<Unit>::AtomIterator& Complex<Unit>::AtomIterator::operator++(int)
{
	operator++();
	return *this;
}

template <class Unit>
typename Complex<Unit>::AtomIterator& Complex<Unit>::AtomIterator::operator=(const typename Complex<Unit>::AtomIterator& rhs)
{
	complex_ = rhs.complex_;
	curr_molecule = rhs.curr_molecule;
	curr_atom = rhs.curr_atom;
	return *this;
}

template <class Unit>
bool Complex<Unit>::AtomIterator::operator==(const typename Complex<Unit>::AtomIterator& rhs)
{
	if(curr_molecule == complex_->molecules_.end())
		return (complex_ == rhs.complex_) && (curr_molecule == rhs.curr_molecule);
	else 
		return (complex_ == rhs.complex_) && (curr_molecule == rhs.curr_molecule) && (curr_atom == rhs.curr_atom);
}

template <class Unit>
bool Complex<Unit>::AtomIterator::operator!=(const typename Complex<Unit>::AtomIterator& rhs)
{
	return !(operator==(rhs));
}

template class Complex<PDBProteinChain>::AtomIterator;

template class Complex<PDBProteinChain>;
template class Complex<Molecule>;

}	//namespace bmpg_uncc_edu::chemistry
}	//namespace bmpg_uncc_edu

#endif

