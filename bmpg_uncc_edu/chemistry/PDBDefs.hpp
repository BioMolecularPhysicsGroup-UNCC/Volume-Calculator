////////////////////////////////////////////////////////////////////////////////
//
// Copyright (c) 2008 by Mike Fairchild
// Author: Mike Fairchild, http://www.mikef.org, mjfairch@gmail.com
//
////////////////////////////////////////////////////////////////////////////////

#ifndef bmpg_uncc_edu_chemistry_PDBDefs_hpp
#define bmpg_uncc_edu_chemistry_PDBDefs_hpp
#include <string>

namespace bmpg_uncc_edu {
namespace chemistry {

using namespace std;

// PDB Files use SPACE when a field is not specified.
// Thus each of the following defaults use a SPACE character.
// Do NOT change these defaults, as these are defined by the
// PDB file standard.
static const char PDB_DEFAULT_REMOTENESS = ' ';
static const char PDB_DEFAULT_ALT_LOCATION = ' ';
static const char PDB_DEFAULT_CHAIN_ID = ' ';
static const char PDB_DEFAULT_BRANCH = ' ';
static const char PDB_DEFAULT_INSERTION_CODE = ' ';

//User defined atoms belong to the following residue and chain to indicate that they are isolated
//static const string NON_RESIDUE_CODE = "---";
//static const string NON_CHAIN_CODE = "-";

} // namespace bmpg_uncc_edu::chemistry
} // namespace bmpg_uncc_edu

#endif

