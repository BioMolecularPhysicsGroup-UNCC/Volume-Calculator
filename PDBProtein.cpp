/**
   Copyright (c) 2008 by Mike Fairchild, Hui Wang
   @author Mike Fairchild, http://www.mikef.org, mjfairch@gmail.com
   @author Hui Wang
*/

#ifndef bmpg_uncc_edu_chemistry_PDBProtein_cpp
#define bmpg_uncc_edu_chemistry_PDBProtein_cpp

#include <sstream>
#include <fstream>
#include <bmpg_uncc_edu/util/Exception.hpp>
//#include <bmpg_uncc_edu/fast/input/EnvironmentCondition.hpp>
#include <bmpg_uncc_edu/chemistry/PDBAtom.hpp>
#include <bmpg_uncc_edu/chemistry/PDBProtein.hpp>
#include <bmpg_uncc_edu/chemistry/PDBProteinResidue.hpp>
#include <bmpg_uncc_edu/chemistry/Complex.hpp>

namespace bmpg_uncc_edu {
namespace chemistry {
using namespace bmpg_uncc_edu::util;

const char all_chains = '-';
	
PDBProtein::PDBProtein() :
		       filename(""),
		       user_id(0),
		       discontiguous_residue_found_(false),
		       decrease_residue_number_found(false)
{
	chain_manager_ = Handle<detail::PDBChainManager>(new detail::PDBChainManager());
	current_chain_ = chain_manager_->new_chain();
}

PDBProtein::~PDBProtein()
{
	
}
	
size_t PDBProtein::num_atoms() const
{
	size_t n = 0;
	residue_iterator_t iter;
	Handle<PDBProteinResidue> r;
	for (iter = residue_begin(); iter != residue_end(); ++iter) {
		r = *iter;
		n += r->num_atoms();
	}
	return n;
}

const PDBAtom * PDBProtein::atom(int n) const
{
	atom_iterator_t iter;
	const PDBAtom * a;
	for (iter = atom_begin(); iter != atom_end(); ++iter) {
		a = *iter;
		if (a->number == n)
			return a;
	}
	return NULL;
}

const PDBAtom * PDBProtein::atom_by_user_id(int id) const
{
	atom_iterator_t iter;
	const PDBAtom * a;
	for (iter = atom_begin(); iter != atom_end(); ++iter) {
		a = *iter;
		if (a->user_id == id)
			return a;
	}
	return NULL;
}

Handle<PDBProteinResidue> PDBProtein::residue(int n,
					      char chain_id,
					      char insertion_code) const
{
	residue_iterator_t iter;
	Handle<PDBProteinResidue> r;
	for (iter = residue_begin(); iter != residue_end(); ++iter) {
		r = *iter;
		// A residue is uniquely specified by BOTH its residue number and its insertion code.
		if (r->number == n && r->chain_id == chain_id && r->insertion_code == insertion_code)
			return r;
	}
	Handle<PDBProteinResidue> null;
	return null;
}

Handle<PDBProteinResidue> PDBProtein::residue_by_user_id(int id) const
{
	residue_iterator_t iter;
	Handle<PDBProteinResidue> r;
	for (iter = residue_begin(); iter != residue_end(); ++iter) {
		r = *iter;
		if (r->number == id)
			return r;
	}
	Handle<PDBProteinResidue> null;
	return null;
}


bool PDBProtein::has_default_chain() const 
{ 
	return (molecules_.find(PDB_DEFAULT_CHAIN_ID) != molecules_.end());
}

bool PDBProtein::has_only_default_chain() const 
{ 
	return (molecules_.size() == 1 && has_default_chain()); 
}

PDBProteinChain * PDBProtein::default_chain() const
{
	if (has_default_chain())
		return molecule(PDB_DEFAULT_CHAIN_ID);
	return NULL;
}

void PDBProtein::build_covalent_bonds()
{
	molecules_map_t::const_iterator iter;
	for (iter = molecules_.begin(); iter != molecules_.end(); ++iter) {
		iter->second->build_covalent_bonds();
	}
}

/**
	Get the covalent bonds of all the molecules in the complex (protein).
*/
PDBProtein::covalent_bonds_t PDBProtein::disulfide_bonds() const
{
	covalent_bonds_t bonds;
	
	molecules_map_t::const_iterator iter;
	for (iter = molecules_.begin(); iter != molecules_.end(); ++iter) {
		const covalent_bonds_t& b = iter->second->disulfide_bonds();
		bonds.insert(bonds.begin(), b.begin(), b.end());
	}
	
	return bonds;
}
	

void PDBProtein::read(const string& fname)
{
	ifstream fin;
	fin.open(fname.c_str());
	if(!fin.is_open()){
		throw util::Exception("PDBProtein::read: failed to open file " + fname, __FILE__,__LINE__);
	}
	
	read(fin);
	fin.close();
	filename = fname;
}

void PDBProtein::read(std::istream & is)
{
	size_t n;
	std::string line;
	
	PDBAtom * atm = NULL;
	Handle<PDBProteinResidue> pres;
	
	while (getline(is,line)) {
		n = line.length();
		 // Line must have at least six characters (PDB record type)
		if (n < 6)
			continue;
		
		// Handle different record types here
		if (line.compare(0, 3, "TER") != 0){// || line.compare(0, 6, "HETATM") == 0){
			atm = read_atom(line);
			crosslink(atm,pres);
		} else {							//new chain was found
			current_chain_ = chain_manager_->new_chain();
		}
		
		// You may implement other record types here later as needed ...
	}
	// Interconnect the file, chains, residues, and atoms together.
}

void PDBProtein::remove_empty_chains()
{
	PDBProtein::molecules_map_t::iterator it = get_empty_chain();
	while(it != molecules_.end()){
		molecules_.erase(it);
	}
}

PDBProtein::molecules_map_t::iterator PDBProtein::get_empty_chain()
{
	
	PDBProtein::molecules_map_t::iterator it;
	for(it = molecules_.begin(); it != molecules_.end(); it++){
		if(it->second->atoms().size() == 0)
			return it;
	}
	return molecules_.end();
}

void PDBProtein::keep_chain_only(const char& id)
{
	PDBProteinChain *chn = molecules_[id];
	
	if(chn == NULL)
		return;
	
	molecules_map_t::iterator it;
	for(it = molecules_.begin(); it != molecules_.end(); it++){
		if(it->first != id){
			delete it->second;
		}
	}
	molecules_.clear();
	molecules_[id] = chn;
}


size_t PDBProtein::num_molecules() const
{
	return molecules_.size();	
}
	
size_t PDBProtein::num_residues() const
{
	size_t total = 0;
	molecules_map_t::const_iterator it;
	for(it = molecules_.begin(); it != molecules_.end(); it++){
		total += it->second->residues.size();
	}
	return total;
}



void PDBProtein::write(std::ostream & os) const
{
	atom_iterator_t aiter;
	const PDBAtom * atm;
	for (aiter = atom_begin(); aiter != atom_end(); ++aiter) {
		atm = *aiter;
		os << atm << endl;
	}
}

PDBAtom* PDBProtein::read_atom(const std::string & line)
{
	// Use auto_ptr for PDBAtom.  Thus if failure or exception occurs, it is cleaned up appropriately.
	
	PDBAtom* atm = NULL;
	
	int ret_val = 0;
	std::auto_ptr<PDBAtom> atom_ptr(new PDBAtom);
	ret_val = atom_ptr->read(line.c_str());
	
	if (ret_val == 0) {
		// Only save this atom record if there were no errors populating it.
		atm = atom_ptr.get();
		atom_ptr.release();						//release auto_ptr ownership of PDBAtom
		return atm;
	}

	return NULL;
}

void PDBProtein::crosslink(PDBAtom* atm, Handle<PDBProteinResidue>& pres)
{
	if(atm == NULL)
		return;
	
	Handle<PDBProteinResidue> res;
	PDBProteinChain * chn = NULL;
	bool new_res_found = false;
	
	atm->chain_id = current_chain_;
	atm->file = this;
	
	chn = molecule(atm->chain_id);
	if (chn == NULL) {
		std::auto_ptr<PDBProteinChain> chn_ptr(new PDBProteinChain);
		chn_ptr->file = this;
		chn_ptr->chain_id = atm->chain_id;
		molecules_[atm->chain_id] = chn_ptr.get();
		chn = chn_ptr.release();
	}
	atm->mol = chn;
	
	res = residue(atm->res_num, atm->chain_id, atm->insertion_code);
	if (res.empty()) {
		Handle<PDBProteinResidue> res_ptr(new PDBProteinResidue);
		res_ptr->number = atm->res_num;
		res_ptr->insertion_code = atm->insertion_code;
		res_ptr->file = this;
		res_ptr->chain = chn;
		res_ptr->chain_id = atm->chain_id;
		res_ptr->aa = atm->aa;
		res = res_ptr;
		chn->residues.push_back(res); 					//assumes all atoms in a residue are on the same chain.
		new_res_found = true;
	}
	res->atoms().push_back(atm);
	chn->atoms().push_back(atm);						//PDBProteinChain (Molecule) owns atoms
	
	atm->res = res;
	
	Handle<PDBProteinResidue> null_res;
	// The following assumes atoms within a residue are listed in increasing 
	// order of atom number.  Same assumption for residues within a chain.
	// I believe this is defined as part of the PDB file format standard.
	if (!pres.empty()) {
		if (pres->number == res->number-1) {
			if(pres->chain_id == res->chain_id){
				pres->next_res = res;
				res->prev_res = pres;
			} else {
				pres->next_res = null_res;
				res->prev_res = null_res;
			}
		}
	}
	if(!pres.empty() && new_res_found && pres->chain_id == res->chain_id && pres->number != res->number - 1 && pres->number >= 0){
		stringstream msg;
		msg << "Jump in residue number in " << pres->number << "-->" << res->number << "\n";
		discontiguous_residue_found_ = true;
	}
	
	if(!pres.empty() && new_res_found && pres->number > res->number){
		decrease_residue_number_found = true;
	}
	
	pres = res;
}

//Handle<EnvironmentCondition>& PDBProtein::set_environment_condition(Handle<EnvironmentCondition> cond)
//{
//	condition_ = cond;
//	condition_->n_repeating_units() = num_residues();
//	return condition_;
//}

std::ostream & operator<<(std::ostream & os, const PDBProtein & file)
{
	file.write(os);
	return os;
}

std::istream & operator>>(std::istream & is, PDBProtein & file)
{
	file.read(is);
	return is;
}



PDBProtein::ResidueIterator::ResidueIterator() : protein_(NULL)
{
}

PDBProtein::ResidueIterator::ResidueIterator(const PDBProtein* p) : protein_(p)
{
	if(protein_->molecules_.size() == 0){
		curr_molecule = protein_->molecules_.end();
		return;
	}
	curr_molecule = protein_->molecules_.begin();
	
	//find a molecule that is neither NULL nor empty
	while (curr_molecule != protein_->molecules_.end() && (curr_molecule->second == NULL || curr_molecule->second->residues.size() == 0)){
			curr_molecule++;
	}

	if(curr_molecule != protein_->molecules_.end()){
		curr_res = curr_molecule->second->residues.begin();
	}

}

PDBProtein::ResidueIterator::ResidueIterator(const PDBProtein* p,
					     molecules_map_t::const_iterator it) :
					     protein_(p),
					     curr_molecule(it)
{
}

PDBProtein::ResidueIterator::ResidueIterator(const PDBProtein* p,
					     molecules_map_t::const_iterator it,
					     residues_t::const_iterator itr) :
					     protein_(p),
					     curr_molecule(it),
					     curr_res(itr)
{
}

PDBProtein::ResidueIterator PDBProtein::residue_begin()
{
	return ResidueIterator(this);
}

PDBProtein::ResidueIterator PDBProtein::residue_end()
{
	return ResidueIterator(this,molecules_.end());
}

PDBProtein::ResidueIterator PDBProtein::residue_begin() const
{
	return ResidueIterator(this);
}

PDBProtein::ResidueIterator PDBProtein::residue_end() const
{
	return ResidueIterator(this,molecules_.end());
}

Handle<PDBProteinResidue> PDBProtein::ResidueIterator::operator*()
{
	return *curr_res;
}

PDBProtein::ResidueIterator& PDBProtein::ResidueIterator::operator++()
{		
	if(curr_res != curr_molecule->second->residues.end()){
		curr_res++;
	}
	if(curr_res == curr_molecule->second->residues.end()){
		curr_molecule++;
		while (curr_molecule != protein_->molecules_.end() && curr_molecule->second->residues.size() == 0){
			curr_molecule++;
		}

		if(curr_molecule != protein_->molecules_.end())
			curr_res = curr_molecule->second->residues.begin();
	}

	return *this;
}

PDBProtein::ResidueIterator& PDBProtein::ResidueIterator::operator++(int)
{
	operator++();
	return *this;
}

PDBProtein::ResidueIterator& PDBProtein::ResidueIterator::operator=(const PDBProtein::ResidueIterator& rhs)
{
	protein_ = rhs.protein_;
	curr_molecule = rhs.curr_molecule;
	curr_res = rhs.curr_res;
	return *this;
}

bool PDBProtein::ResidueIterator::operator==(const PDBProtein::ResidueIterator& rhs)
{
	if(curr_molecule == protein_->molecules_.end())
		return (protein_ == rhs.protein_) && (curr_molecule == rhs.curr_molecule);
	else 
		return (protein_ == rhs.protein_) && (curr_molecule == rhs.curr_molecule) && (curr_res == rhs.curr_res);
}

bool PDBProtein::ResidueIterator::operator!=(const PDBProtein::ResidueIterator& rhs)
{
	return !(operator==(rhs));
}


} // namespace bmpg_uncc_edu::chemistry
} // namespace bmpg_uncc_edu

#endif

