/**
  Copyright (c) 2008 by Mike Fairchild, Hui Wang
  @author: Mike Fairchild, http://www.mikef.org, mjfairch@gmail.com
  @author Hui Wang
*/

#ifndef bmpg_uncc_edu_chemistry_PDBAtom_cpp
#define bmpg_uncc_edu_chemistry_PDBAtom_cpp

#include <string.h>
#include <cstdlib>
#include <iostream>
#include <sstream>
#include <iomanip>
#include <bmpg_uncc_edu/util/Exception.hpp>
#include <bmpg_uncc_edu/util/logger/LoggerFactory.hpp>
#include <bmpg_uncc_edu/chemistry/PDBDefs.hpp>
#include <bmpg_uncc_edu/chemistry/Element.hpp>
#include <bmpg_uncc_edu/chemistry/PDBAtom.hpp>
#include <bmpg_uncc_edu/chemistry/PDBProteinChain.hpp>
#include <bmpg_uncc_edu/chemistry/PDBProtein.hpp>
#include <bmpg_uncc_edu/chemistry/helper/PDBHelper.hpp>

namespace bmpg_uncc_edu {
namespace chemistry {
using namespace std;
using namespace bmpg_uncc_edu::chemistry::helper;
using namespace bmpg_uncc_edu::util::logger;

PDBAtom::PDBAtom()
{
	clear();
}

PDBAtom::~PDBAtom()
{
	clear();
}

void PDBAtom::clear()
{
	// Clear defined members
	number = 0;
	atom_name[0] = atomic_symbol[0] = res_sname[0] = record_id[0] = segment_id[0] = charge[0] = '\0';
	remoteness = PDB_DEFAULT_REMOTENESS;
	branch = PDB_DEFAULT_BRANCH;
	alt_location = PDB_DEFAULT_ALT_LOCATION;
	chain_id = PDB_DEFAULT_CHAIN_ID;
	insertion_code = PDB_DEFAULT_INSERTION_CODE;
	temp_factor = occupancy = 0;
	x = y = z = 0;
	res_num = 0;
	user_id = 0;
	aa = NULL;
	file = NULL;
	mol = NULL;
}

bool PDBAtom::on_default_chain() const
{
	return (chain_id == PDB_DEFAULT_CHAIN_ID);
}

bool PDBAtom::on_branch() const
{
	return (branch != PDB_DEFAULT_BRANCH);
}

bool PDBAtom::is_heavy() const
{
	return (strncmp(atomic_symbol, " H", 2) != 0);
}

int PDBAtom::read(std::istream & is)
{
	if (!is || is.eof())
		return ERR_EOF;
	
	// Read a line from the stream (up to 127 characters + NUL).
	char str[128];
	is.getline(str, 128);
	
	return read(str);
}

int PDBAtom::read(const char * str)
{
	size_t n;
	char tmp[32];
	char buf[128];
	
	AminoAcidLibrary * aa_db = AminoAcidLibrary::instance();
	
	// Setup buffer pointer.  It points to the user-supplied string, provided
	// that string has at least 80 characters.  If it doesn't, then copy
	// contents of user-supplied string and pad with spaces to 80 characters.
	const char * buffer = str;
	n = strlen(str);
	if (n < 80) {
		strncpy(buf, str, n);
		memset(&buf[n], ' ', 80-n);
		buf[80]='\0';
		buffer = buf;
	}
	
	if (strncmp(buffer, "TER", 3) == 0)					//TER in record
		return TER_RECORD_FOUND;
	
	
	// A PDB ATOM record starts with "ATOM"
	if (strncmp(buffer, "ATOM", 4) != 0)// && strncmp(buffer, "HETATM", 6) != 0)
		return ERR_NOT_ATOM;
	
	// Defined member
	strncpy(tmp, &buffer[6], 5); tmp[5]='\0'; number = atoi(tmp);
	strncpy(atom_name, &buffer[12], 4); atom_name[4] = '\0';

	if (strncmp(atom_name, " UNK", 4) == 0)
		return ERR_UNKNOWN_ATOM;

	element_ = PDBHelper::atom_to_element(atom_name);
	
	alt_location = buffer[16];
	if(alt_location != 'A' && alt_location != PDB_DEFAULT_ALT_LOCATION)	//skip atoms whoes alternative location is not A or ' '
		return ERR_NOT_DEFAULT_ALT_LOCATION;
	
	strncpy(res_sname, &buffer[17], 3); res_sname[3] = '\0';
	chain_id = buffer[21];
	strncpy(tmp, &buffer[22], 4); tmp[4]='\0'; res_num = atoi(tmp);
	insertion_code=buffer[26];
	strncpy(tmp, &buffer[30], 8); tmp[8]='\0'; x = atof(tmp);
	strncpy(tmp, &buffer[38], 8); tmp[8]='\0'; y = atof(tmp);
	strncpy(tmp, &buffer[46], 8); tmp[8]='\0'; z = atof(tmp);
	strncpy(tmp, &buffer[54], 6); tmp[6]='\0'; occupancy = atof(tmp);
	strncpy(tmp, &buffer[60], 6); tmp[6]='\0'; b_value = atof(tmp);
	strncpy(record_id, &buffer[72], 8); record_id[8]='\0';
	strncpy(segment_id, &buffer[72], 4); segment_id[4]='\0';
	strncpy(charge, &buffer[78], 2); charge[2]='\0';

	// Derived members
	strncpy(atomic_symbol, &buffer[12], 2); atomic_symbol[2] = '\0';
	remoteness = buffer[14];
	branch = buffer[15];

	aa = aa_db->find_by_sname(res_sname);
	
	if(aa == NULL){
		Logger* logger = LoggerFactory::default_logger();
		stringstream msg;
		msg << "Residue " << res_sname << " is not found in amino acid library;" << " Atom number : " << number << ".";
		logger->info(msg.str());
		return ERR_INVALID_RESIDUE;
	}
	
	return 0;
}

void PDBAtom::write(std::ostream & os, bool verbose) const
{
	using namespace std;
	
	if (!verbose) {
		os << setfill(' ');
		os << setw(6) << left << "ATOM";
		os << setw(5) << right << number << ' ';
		os << setw(4) << atom_name;
		os << alt_location;
		os << setw(3) << res_sname << ' ';
		os << chain_id;
		os << setw(4) << res_num;
		os << insertion_code;
		os << "   ";
		os << setw(8) << fixed << setprecision(3) << x;
		os << setw(8) << setprecision(3) << y;
		os << setw(8) << setprecision(3) << z;
		os << setw(6) << setprecision(2) << occupancy;
		os << setw(6) << b_value;
		os << "      ";
		os << setw(8) << record_id;
	}
	else {
		// Verbose mode here
		os << "User ID=" << user_id << endl;
		os << "Number=" << number << endl;
		os << "Atom Name=\"" << atom_name << "\"" << endl;
		os << "Alt Location='" << alt_location << "'" << endl;
		os << "Residue Short Name=\"" << res_sname << "\"" << endl;
		os << "Chain ID='" << chain_id << "'" << endl;
		os << "Residue number=" << res_num << endl;
		os << "Insertion code='" << insertion_code << "'" << endl;
		os << "X=" << x << endl;
		os << "Y=" << y << endl;
		os << "Z=" << z << endl;
		os << "Occupancy=" << occupancy << endl;
		os << "Temp factor/B value=" << temp_factor << endl;
		os << "Record ID=\"" << record_id << "\"" << endl;
		os << "Segment ID=\"" << segment_id << "\"" << endl;
		os << "Charge=\"" << charge << "\"" << endl;
		os << "Atomic Symbol=\"" << atomic_symbol << "\"" << endl;
		os << "Remoteness='" << remoteness << "'" << endl;
		os << "Branch='" << branch << "'" << endl;
		os << "On default chain? " << (on_default_chain() ? "Yes" : "No") << endl;
		os << "On a branch? " << (on_branch() ? "Yes" : "No") << endl;
		os << "Is heavy? " << (is_heavy() ? "Yes" : "No") << endl;
	}
	cout << fixed << setprecision(16);
}


PDBAtom::PDBAtom(const PDBAtom& rhs) :Atom(rhs)
{
	operator=(rhs);
}

PDBAtom& PDBAtom::operator=(const PDBAtom& rhs)
{
	number = rhs.number;
	strncpy(atom_name,rhs.atom_name,5);
	alt_location = rhs.alt_location;
	strncpy(res_sname,rhs.res_sname,4);
	chain_id = rhs.chain_id;
	res_num = rhs.res_num;
	insertion_code = rhs.insertion_code;
	occupancy = rhs.occupancy;
	
	temp_factor = rhs.temp_factor;
	b_value = rhs.b_value;
	strncpy(record_id,rhs.record_id,9);
	strncpy(segment_id,rhs.segment_id,5);
	strncpy(charge,rhs.charge,3);
	strncpy(atomic_symbol,rhs.atomic_symbol,3);
	
	remoteness = rhs.remoteness;
	
	branch = rhs.branch;
	x = rhs.x;
	y = rhs.y;
	z = rhs.z;
	user_id = rhs.user_id;

	file = rhs.file;
	mol = rhs.mol;
	res = rhs.res;
	
	aa = rhs.aa;
	element_ = rhs.element_;
	return *this;
}

std::ostream & operator<<(std::ostream & os, const PDBAtom & atom)
{
	atom.write(os, false);
	return os;
}

std::istream & operator>>(std::istream & is, PDBAtom & atom)
{
	atom.read(is);
	return is;
}

}	//namespace bmpg_uncc_edu::chemistry
}	//namespace bmpg_uncc_edu

#endif



