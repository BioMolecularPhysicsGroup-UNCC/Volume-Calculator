/**
   Copyright (c) 2008 by Hui Wang
   @author Hui Wang
   @date Jan 20, 2009
*/


#ifndef bmpg_uncc_edu_chemistry_CovalentBonder_cpp
#define bmpg_uncc_edu_chemistry_CovalentBonder_cpp

#include <iostream>
#include <sstream>
#include <vector>
#include <cmath>
#include <bmpg_uncc_edu/util/Exception.hpp>
#include <bmpg_uncc_edu/chemistry/library/BondDistanceTable.hpp>
#include <bmpg_uncc_edu/chemistry/CovalentBonder.hpp>
#include <bmpg_uncc_edu/chemistry/CovalentBond.hpp>
#include <bmpg_uncc_edu/chemistry/library/BondDistanceTable.hpp>
#include <bmpg_uncc_edu/chemistry/PDBAtom.hpp>
#include <bmpg_uncc_edu/algorithms/HashGrid.hpp>
#include <bmpg_uncc_edu/chemistry/PDBProtein.hpp>
#include <bmpg_uncc_edu/chemistry/Molecule.hpp>
#include <bmpg_uncc_edu/chemistry/helper/PDBHelper.hpp>
#include <bmpg_uncc_edu/util/StringUtil.hpp>
#include <bmpg_uncc_edu/util/logger/LoggerFactory.hpp>

namespace bmpg_uncc_edu {
namespace chemistry {
using namespace bmpg_uncc_edu::algorithms;
using namespace bmpg_uncc_edu::chemistry::helper;
using namespace bmpg_uncc_edu::util;
using namespace bmpg_uncc_edu::util::logger;

double CovalentBonder::hard_distance_cutoff_ = 6.0;

bool CovalentBonder::meets_atom_atom_distance_cutoff(const PDBAtom* atom1,
						     const PDBAtom* atom2)
{  
	const Element *elem1 = atom1->element();
	const Element *elem2 = atom2->element();

	BondDistanceTable * table = BondDistanceTable::instance();
	double cutoff = table->find(elem1, elem2);
	return (PDBHelper::distance(atom1, atom2) <= cutoff );
}


CovalentBonder::cov_bonds_t CovalentBonder::build_covalent_bonds(PDBProtein& protein)
{
	disulfide_bonds_.clear();
	std::vector<PDBProteinChain*> chains = protein.molecules();

	CovalentBonder::cov_bonds_t bonds;
	for(size_t i = 0; i < chains.size(); i++){
		CovalentBonder::cov_bonds_t bs = build_covalent_bonds(*chains.at(i));
		bonds.insert(bonds.begin(),bs.begin(),bs.end());
	}
	built_ = true;
	return bonds;
}


CovalentBonder::cov_bonds_t CovalentBonder::build_covalent_bonds(Molecule & mol)
{
	Molecule::atoms_t & atoms = mol.atoms();
	
	HashGrid grid(atoms, grid_length_);					//set up the grid for hashing
	cov_bonds_t bonds;
	bool passes_distance_check, templates_confirm_neighbors, templates_confirm_nonneighbors, cannot_discern_bond_from_templates;
	
	Logger * logger = LoggerFactory::default_logger();
	
	Molecule::atoms_t::iterator it,it2;
	
	for(it = atoms.begin(); it != atoms.end(); it++){			//for each atom in the molecule
		const PDBAtom *atom1 = *it;
		Molecule::atoms_t natoms = grid.atom_neighbors(atom1);		//find all the possible neighbors of the atom
		for(it2 = natoms.begin(); it2 != natoms.end(); it2++){		//for each neighbor
			const PDBAtom *atom2 = *it2;

			if(atom2->number <= atom1->number)continue;

			passes_distance_check = meets_atom_atom_distance_cutoff(atom1,atom2);
			string s1 = atom1->atom_name;
			string s2 = atom2->atom_name;
						 
			
			if(atom1->res == atom2->res){				//same residue
				const AminoAcid *aa = atom1->aa;
				templates_confirm_neighbors = aa->has_atom(atom1->atom_name) && aa->has_atom(atom2->atom_name) &&
							      aa->is_neighbor(atom1->atom_name, atom2->atom_name);
				
				templates_confirm_nonneighbors = aa->has_atom(atom1->atom_name) && aa->has_atom(atom2->atom_name) &&
							         !templates_confirm_neighbors;
				cannot_discern_bond_from_templates = !templates_confirm_neighbors && !templates_confirm_nonneighbors;
				
				if(templates_confirm_neighbors){
					if(PDBHelper::distance(atom1,atom2) > hard_distance_cutoff_){
						stringstream msg;
						msg << "The covalent bond between " << *atom1 << " and " << *atom2 << " in residue " << aa->short_name() << 
								" is defined in the amino acid library but exceeds the hard distance cutoff, d = " <<
							PDBHelper::distance(atom1, atom2) << ", cutoff = " << hard_distance_cutoff_ << " Angstroms.";
						logger->error(msg.str());
					}
					logger->debug("template bond 1");
					
					bonds.push_back(new CovalentBond(atom1,atom2));
				} else if(cannot_discern_bond_from_templates && passes_distance_check){
					logger->debug("distance bond 2");
					bonds.push_back(new CovalentBond(atom1,atom2));
				}
				 
				if(templates_confirm_neighbors && !passes_distance_check){
					stringstream msg;
					msg << "The covalent bond between " << *atom1 << " and " << *atom2 << " in residue " << aa->short_name() << 
						 		" is defined in the amino acid library but exceeds the bond distance cutoff, d = "
								<< PDBHelper::distance(atom1, atom2) << " Angstroms";
					logger->warn(msg.str());
					bonds.push_back(new CovalentBond(atom1,atom2));
				}
				
			}	//end intra-residue check
			else {	//inter-residue
				if(passes_distance_check){
					if(check_peptide_bond(atom1,atom2)){					//peptide bond
						logger->debug("peptide bond 3");
						bonds.push_back(new CovalentBond(atom1,atom2));
					} else if (check_disulfide_bond(atom1,atom2)){				//disulfide bond
						logger->debug("disulfide bond 4");
						Handle<CovalentBond> bond = new CovalentBond(atom1,atom2);
						bonds.push_back(bond);
						disulfide_bonds_.push_back(bond);
					} else if (check_internucleic_bond(atom1,atom2)){			//inter-nuclei bond
						logger->debug("internucleic bond 5");
						bonds.push_back(new CovalentBond(atom1,atom2));
					} else {
						stringstream msg;
						msg << "Atoms " << *atom1 << " and " << *atom2 << " in residue " << atom1->aa->short_name() << 
						 		" is within bond distance but not recognized as peptide/disulfide/inter-nucleic bonds.";
						logger->warn(msg.str());
					}
					
				}	//end of passes_distance_check
				
			}		//end of inter-residue check
		}	//end of looping over neighbors
		
	}		//end of looping over atoms
	
	built_ = true;
	return bonds;
}

CovalentBonder::disulfide_bonds_t CovalentBonder::disulfide_bonds() const
{
	if(!built_){
		throw Exception("No covalent bonds have been built yet for any molecule.", __FILE__, __LINE__);
	}
	return disulfide_bonds_;
}
	

bool CovalentBonder::check_peptide_bond(const PDBAtom * atom1,
					const PDBAtom * atom2)
{
	string s1 = atom1->atom_name;
	string s2 = atom2->atom_name;
	return (s1 == " N  " && s2 == " C  ") || (s2 == " N  " && s1 == " C  ");
}


bool CovalentBonder::check_disulfide_bond(const PDBAtom * atom1,
					  const PDBAtom * atom2)
{
	string s1 = atom1->atom_name;
	string s2 = atom2->atom_name;
	string res1 = atom1->aa->short_name();
	string res2 = atom2->aa->short_name();
	
	return (s1 == " SG " && s2 == " SG ") && (res1 == "CYS" || res1 == "CYX") && (res2 == "CYS" || res2 == "CYX");
}

bool CovalentBonder::check_internucleic_bond(const PDBAtom * atom1,
					     const PDBAtom * atom2)
{
	string s1 = atom1->atom_name;
	string s2 = atom2->atom_name;
	return (s1 == " O3 " && s2 == " P  ") || (s2 == " O3 " && s1 == " P  ");
}

}	//namespace bmpg_uncc_edu::chemistry
}	//namespace bmpg_uncc_edu

#endif

