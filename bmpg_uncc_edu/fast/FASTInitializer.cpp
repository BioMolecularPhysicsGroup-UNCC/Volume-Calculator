/**
   Copyright (c) 2008 by Hui Wang
   @author Hui Wang
   Modified by Jenny Farmer 2019
   @date Jan 20, 2009
*/


#ifndef bmpg_uncc_edu_fast_FASTInitializer_cpp
#define bmpg_uncc_edu_fast_FASTInitializer_cpp

#include "FASTInitializer.hpp"

#include <fstream>
#include <bmpg_uncc_edu/chemistry/library/AminoAcidLibrary.hpp>
#include <bmpg_uncc_edu/chemistry/library/PeriodicTable.hpp>

using namespace std;
using namespace bmpg_uncc_edu::chemistry::library;


FASTInitializer::FASTInitializer() {
}

void FASTInitializer::load_user_defined_libraries()
{	
	string  pt_filename = "periodic_table.csv";
        string  aa_db_filename = "amino_acid_term.xml";
	string  path_to_metadata = "";	
        
        ifstream fin;
        
        //load periodic table
	PeriodicTable *pd = PeriodicTable::instance();
	pt_filename = path_to_metadata + pt_filename;	
	fin.open(pt_filename.c_str());
	pd->load(fin);
	fin.close();
	
	//load amino acid library
	aa_db_filename = path_to_metadata + aa_db_filename;
	AminoAcidLibrary * aa_db = AminoAcidLibrary::instance();
	aa_db->load(aa_db_filename.c_str());
	
}

void FASTInitializer::close() {
}


#endif
