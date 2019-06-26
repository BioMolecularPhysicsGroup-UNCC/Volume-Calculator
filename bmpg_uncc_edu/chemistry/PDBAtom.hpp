////////////////////////////////////////////////////////////////////////////////
//
// Copyright (c) 2008 by Mike Fairchild
// Author: Mike Fairchild, http://www.mikef.org, mjfairch@gmail.com
//
////////////////////////////////////////////////////////////////////////////////

#ifndef bmpg_uncc_edu_pdb_PDBAtom_hpp
#define bmpg_uncc_edu_pdb_PDBAtom_hpp

#include <string>
#include <bmpg_uncc_edu/util/Handle.hpp>
#include <bmpg_uncc_edu/chemistry/Atom.hpp>
#include <bmpg_uncc_edu/chemistry/library/AminoAcidLibrary.hpp>

// See (OFFICIAL): http://tinyurl.com/5souud
// See (UNOFFICIAL): http://tinyurl.com/5mmb42

// Note: Although the new PDB standard requires columns 77-78 of an ATOM record
// to contain the two-character element symbol, many existing PDB files do
// not have this information here.  Therefore, it is safest to extract these
// two characters from the first two characters of the atom field, which is
// (almost) always the two character element symbol.  The main exception to
// this happens in hydrogens, where the first character is sometimes a digit.

// Why am I not using "good OO design" by hiding data members of the PDBAtom
// class and providing get/set accessor methods?  Answer: This is good OO
// design.  PDBAtom, PDBResidue, and PDBChain are basically POD (Plain Old
// Data). For POD (like structs), it's best to let the user have direct access.
// Hence the public members.  I still make protected/private implementation
// details.

namespace bmpg_uncc_edu {
namespace chemistry {
using namespace bmpg_uncc_edu::util;
using namespace bmpg_uncc_edu::chemistry::library;

class PDBProtein;								// Forward declaration
class PDBProteinChain;								// Forward declaration
class PDBProteinResidue;							// Forward declaration
class Molecule;
class Element;


/**
 * PDBAtom represents the information found in a PDB file for a single atom. 
 */

class PDBAtom : public Atom
{
public:
	// Error codes, signifying problems when populating a PDBAtom object.
	static const int ERR_EOF = -1;
	static const int ERR_NOT_ATOM = -2;
	static const int ERR_INVALID_REMOTENESS = -3;
	static const int ERR_INVALID_RESIDUE = -4;

	static const int ERR_NOT_DEFAULT_ALT_LOCATION = -5;
	static const int ERR_UNKNOWN_ATOM = -6;
	static const int TER_RECORD_FOUND = -7;
	
	// Defined members of a PDB ATOM record, according to the PDB standard
	int number;								// atom number in the PDB file
	char atom_name[4+1];							// atom_name[0-1]=element symbol, atom_name[2]=remoteness, atom_name[3]=branch)
	char alt_location;
	char res_sname[3+1];							// 3 character short name for residue
	char chain_id;								// together with res_num and insertion_code uniquely identifies a residue in a PDB file
	int res_num;								// together with insertion_code uniquely identifies a residue on a chain
	char insertion_code;							// together with res_num uniquely identifies a residue on a chain
	float occupancy;
	union {
		float temp_factor;						// The same location is used for either the
		float b_value;							// temp factor or the B value.  Hence the union.
	};
	char record_id[8+1];
	char segment_id[4+1];
	char charge[2+1];
	
	// Derived members (i.e. derived from the defined members listed above)
	char atomic_symbol[2+1];						// derived from atom member
	char remoteness;							// derived from atom member
	char branch;								// derived from atom member
		
	// Miscellanoues members
	int user_id;								// User-defined atom id.  Define as you please
	const AminoAcid * aa;							// derived from res_sname member
	PDBProtein * file;							// Parent PDBFile to which this atom belongs
	Molecule* mol;								// Parent PDBChain to which this atom belongs
	Handle<PDBProteinResidue> res;						// Parent PDBResidue to which this atom belongs

	const Element* atom_to_element(const char* aname) const;

	// Constructors, destructors, memory mangement
	PDBAtom();
	~PDBAtom();
	void clear();
	
	// Query methods
	bool on_default_chain() const;						// True if not on the default chain
	bool on_branch() const;							// True if on a branch (e.g. ring or tree)
	bool is_heavy() const;							// True (by definition) if a hydrogen atom

	// Methods for populating this object
	int read(std::istream & is);
	int read(const char * str);
	int read(const std::string & str) { return read(str.c_str()); }
	void write(std::ostream & ostr, bool verbose = false) const;
	
	PDBAtom & operator=(const PDBAtom & rhs);				// hide assignment
	PDBAtom(const PDBAtom & rhs);						// hide copy constructor

	
	// Friend functions
	friend std::ostream & operator<<(std::ostream & os, const PDBAtom & atom);
	friend std::istream & operator>>(std::istream & is, PDBAtom & atom);

	
};

}	//bmpg_uncc_edu::chemistry
}	//bmpg_uncc_edu

#endif

