/**
   Copyright (c) 2008 by Mike Fairchild, Hui Wang
   @author: Mike Fairchild, http://www.mikef.org, mjfairch@gmail.com
   @author: Hui Wang
*/


#ifndef bmpg_uncc_edu_chemistry_PDBProteinResidue_cpp
#define bmpg_uncc_edu_chemistry_PDBProteinResidue_cpp

#include <string.h>
#include <cstdlib>
#include <boost/functional/hash.hpp>
#include <bmpg_uncc_edu/fast/ParameterFile.hpp>
#include <bmpg_uncc_edu/chemistry/PDBAtom.hpp>
#include <bmpg_uncc_edu/chemistry/helper/PDBHelper.hpp>
#include <bmpg_uncc_edu/chemistry/PDBProteinResidue.hpp>
#include <bmpg_uncc_edu/chemistry/PDBProtein.hpp>
#include <bmpg_uncc_edu/chemistry/library/AminoAcidLibrary.hpp>

namespace bmpg_uncc_edu {
namespace chemistry {
using namespace std;
using namespace bmpg_uncc_edu::chemistry::helper;

PDBProteinResidue::PDBProteinResidue() :
				     aa(NULL),
				     chain(NULL),
				     file(NULL),
				     chain_id(PDB_DEFAULT_CHAIN_ID),
				     number(0), 
				     insertion_code(PDB_DEFAULT_INSERTION_CODE),
				     user_id(0),
				     inter_cvdw_energy(0),
				     inter_hbond_energy(0),
				     nnc_(0),
				     asa_(0),
				     sec_struct_("NULL"),
				     total_native_energy_(0)
{
	// Reserve space for an "average" number of atoms in a PDB residue.
	// Admittedly this is a guess.  The data structure will grow
	// automatically as needed.  This call to reserve() improves performance
	// at the beginning.
	atoms_.reserve(10);
}

PDBProteinResidue::~PDBProteinResidue()
{
	atoms_.clear();
}
 

PDBAtom * PDBProteinResidue::atom(const char * name,
				  bool try_alt_hydrogen_names) const
{
	if (strlen(name) < 4)
		return NULL;
	
	double b = 0;
	PDBAtom* a = NULL;
	atoms_t as = atoms_with_all_occupancy(name, try_alt_hydrogen_names);
	
	if(atoms_.size() == 0){
		return NULL;
	} else if(as.size() == 1){		
		return as.at(0);
	} else {
		bmpg_uncc_edu::fast::ParameterFile* param_file = bmpg_uncc_edu::fast::ParameterFile::instance();
		if(!param_file->atom_by_occupancy()){
			if(as.size() > 0)
				return as.at(0);
		} else {
			atoms_t::const_iterator iter;
			for(iter = as.begin(); iter != as.end(); iter++){
				if( (*iter)->alt_location == 'A' ){					
					return (*iter);
				}
				if( (*iter)->occupancy > b){
					a = *iter;
					b = (*iter)->occupancy;
				}
				
			}
		}
	}
	return a;
}

bool PDBProteinResidue::has_atom(const std::string& s,
				 bool try_alt_hydrogen_names) const
{
	return (atom(s,try_alt_hydrogen_names) != NULL);	
}
	
                        
PDBProteinResidue::atoms_t PDBProteinResidue::atoms_with_all_occupancy(const char * name,
								       bool try_alt_hydrogen_names) const
{
	atoms_t result;
	
	atoms_t::const_iterator iter;
	PDBAtom * a = NULL;
	const char * alternate_name = NULL;
	for (iter = atoms_.begin(); iter != atoms_.end(); ++iter) {
		a = *iter;
		if (strncmp(name, a->atom_name, 4) == 0)
			result.push_back(a);
	}

	// If we're here the lookup failed.  If the user named a hydrogen atom and
	// wants to try alternate hydrogen naming rules, then do so.
	if (name[1] == 'H' && try_alt_hydrogen_names) {
		if ((alternate_name = PDBHelper::alternate_hydrogen_name(name)) != NULL) {
			for (iter = atoms_.begin(); iter != atoms_.end(); ++iter) {
				a = *iter;
				if (strncmp(alternate_name, a->atom_name, 4) == 0)
					result.push_back(a);
			}
		}
	}

	return result;
}


PDBAtom * PDBProteinResidue::atom(const char * symbol,
				  char remoteness,
				  char branch,
				  bool try_alt_hydrogen_names) const
{
	if (strlen(symbol) < 2)
		return NULL;
	char name[5];
	name[0] = symbol[0];
	name[1] = symbol[1];
	name[2] = remoteness;
	name[3] = branch;
	name[4] = '\0';
	return atom(name, try_alt_hydrogen_names);
}


size_t PDBProteinResidue::hash_code()
{
	boost::hash<PDBProteinResidue> key_hasher;
	return key_hasher(*this);	
}

size_t hash_value(const PDBProteinResidue& res)
{
	boost::hash<string> hasher;
	string s =  PDBHelper::residue_info(res);
	return hasher(s);	
}

std::ostream & operator<<(std::ostream & os, const PDBProteinResidue & res)
{
	PDBProteinResidue::atoms_t::const_iterator aiter;
	const PDBAtom * atm;
	for (aiter = res.atoms_.begin(); aiter != res.atoms_.end(); ++aiter) {
		atm = *aiter;
		os << *atm;
	}
	return os;
}
	

}	//namespace bmpg_uncc_edu::chemistry
}	//namespace bmpg_uncc_edu

#endif

