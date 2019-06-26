/**
   Copyright (c) 2008 by Hui Wang, Chuanbin Du
   @author Hui Wang, Chuanbin Du
   @date Jun 03, 2009
*/


#ifndef bmpg_uncc_edu_chemistry_library_loader_AALibLoader_cpp
#define bmpg_uncc_edu_chemistry_library_loader_AALibLoader_cpp

#include <cmath>
#include <fstream>
#include <sstream>
#include <bmpg_uncc_edu/util/logger/LoggerFactory.hpp>
#include <bmpg_uncc_edu/chemistry/library/AminoAcidLibrary.hpp>
#include <bmpg_uncc_edu/chemistry/library/loader/AALibLoader.hpp>
#include <bmpg_uncc_edu/util/Exception.hpp>
#include <bmpg_uncc_edu/util/xml/xmlParser.hpp>

namespace bmpg_uncc_edu {
namespace chemistry {
namespace library {
namespace loader {
using namespace std;
using namespace bmpg_uncc_edu::util;
using namespace bmpg_uncc_edu::util::xml;
using namespace bmpg_uncc_edu::util::logger;

AALibLoader::AALibLoader(const string& fname,
			 AminoAcidLibrary& lib) :
			 fname_(fname),
			 lib_(lib)
{
	
}

void AALibLoader::load()
{	
	XMLNode xMainNode;
	XMLNode xNode_AminoAcid_i;
	XMLNode xNode_AminoAcid_i_Covalent;
	XMLNode xNode_AminoAcid_i_Covalent_j;
			
	xMainNode = XMLNode::openFileHelper(fname_.c_str());
	
	XMLNode xNode = xMainNode.getChildNode("AminoAcidLibrary");

	int nAminoAcid = xNode.nChildNode("AminoAcid");
	
	string aa_long_name, aa_formula;
	string atom1, atom2;
	for(int i = 0; i < nAminoAcid; i++){		
		xNode_AminoAcid_i = xNode.getChildNode("AminoAcid",i);
		string aa_short_name = xNode_AminoAcid_i.getAttribute("short_name");
		
		AminoAcid* aa = new AminoAcid();
		
		//all upper cases for amino acid names
		std::transform(aa_short_name.begin(), aa_short_name.end(), aa_short_name.begin(), ::toupper);
	
		aa->short_name() = aa_short_name;
		aa_long_name = xNode_AminoAcid_i.getAttribute("long_name");
		std::transform(aa_long_name.begin(), aa_long_name.end(), aa_long_name.begin(), ::toupper);

		aa->long_name() = aa_long_name;
		aa_formula = xNode_AminoAcid_i.getAttribute("formula");
		aa->formula() = aa_formula;
		
		string fed_calculation(xNode_AminoAcid_i.getAttribute("fed_calculation"));
		
		
		string s_acid_base(xNode_AminoAcid_i.getAttribute("acid_base"));
		if(s_acid_base.compare("neutral") == 0){
			aa->acid_base = AminoAcid::neutral;
		} else if(s_acid_base.compare("acid") == 0){
			aa->acid_base = AminoAcid::acid;
		} else if(s_acid_base.compare("base") == 0){
			aa->acid_base = AminoAcid::base;
		} else {
			throw Exception("Unrecognized acid base type for amino acid " + aa->short_name(), __FILE__, __LINE__);
		}
		
		bool aa_for_fed = fed_calculation.compare("yes") == 0;
		string s_aa_include(xNode_AminoAcid_i.getAttribute("include"));
		bool aa_include = s_aa_include.compare("yes") == 0;
		
		if(lib_.has(aa_short_name) && aa_include){
			Logger* logger = LoggerFactory::default_logger();
			logger->critical("Amino acid " + aa_short_name + " is already in the amino acid library.");
		}
		
		
		XMLNode node_solvation = xNode_AminoAcid_i.getChildNode("SolvationProperty");
		string s_surface = node_solvation.getAttribute("flanking_asa");
		stringstream ssd(s_surface);
		double flanking_asa;
		ssd >> flanking_asa;
		aa->flanking_asa() = flanking_asa;
		
		s_surface = node_solvation.getAttribute("naked_asa");
		stringstream ssn(s_surface);
		double naked_asa;
		ssn >> naked_asa;
		aa->naked_asa() = naked_asa;
		
		
		xNode_AminoAcid_i_Covalent = xNode_AminoAcid_i.getChildNode("CovalentConnections");		
		int nConnection = xNode_AminoAcid_i_Covalent.nChildNode("Connection");

	
		for(int j = 0; j < nConnection; j++){
			xNode_AminoAcid_i_Covalent_j = xNode_AminoAcid_i_Covalent.getChildNode("Connection", j);
			atom1 = xNode_AminoAcid_i_Covalent_j.getAttribute("atom1");
			atom2 = xNode_AminoAcid_i_Covalent_j.getAttribute("atom2");
			aa->add_connection(atom1,atom2);	
		}
		
		//Read in hybridization states
		XMLNode node_protonation = xNode_AminoAcid_i.getChildNode("Protonation");		
		int n_prot_atoms = node_protonation.nChildNode("Protonation_atom");
		for(int j = 0; j < n_prot_atoms; j++){
			XMLNode node_prot_atom = node_protonation.getChildNode("Protonation_atom", j);
			string atom_name = node_prot_atom.getAttribute("name");
			
			int n_hyb_states = node_prot_atom.nChildNode("Hybridization");
			for(int k = 0; k < n_hyb_states; k++){
				XMLNode node_hyb = node_prot_atom.getChildNode("Hybridization",k);
				string hybridization = node_hyb.getAttribute("value");
				aa->add_hybridization(atom_name, hybridization);
			}
		}
		
		if(aa_include)
			lib_.add(aa);

		if(aa_for_fed){
			lib_.add_amino_acid_to_fed(aa->short_name());
		}
		
	}
}

}	//namespace bmpg_uncc_edu::chemistry::library::loader	
}	//namespace bmpg_uncc_edu::chemistry::library	
}	//namespace bmpg_uncc_edu::chemistry
}	//namespace bmpg_uncc_edu

#endif

