/**
   Copyright (c) 2008 by Hui Wang
   @author Hui Wang
   @date Jan 20, 2009
*/


#ifndef bmpg_uncc_edu_chemistry_library_AminoAcidLibrary_cpp
#define bmpg_uncc_edu_chemistry_library_AminoAcidLibrary_cpp

#include <iomanip>
#include <cmath>
#include <bmpg_uncc_edu/util/logger/LoggerFactory.hpp>
#include <bmpg_uncc_edu/chemistry/Element.hpp>
#include <bmpg_uncc_edu/chemistry/library/AminoAcidLibrary.hpp>
#include <bmpg_uncc_edu/chemistry/library/loader/AALibLoader.hpp>
#include <bmpg_uncc_edu/chemistry/library/PeriodicTable.hpp>
#include <bmpg_uncc_edu/chemistry/helper/PDBHelper.hpp>
#include <bmpg_uncc_edu/util/Exception.hpp>

namespace bmpg_uncc_edu {
namespace chemistry {
namespace library {
using namespace std;
using namespace bmpg_uncc_edu::util;
using namespace bmpg_uncc_edu::chemistry::helper;
using namespace bmpg_uncc_edu::chemistry::library::loader;
using namespace bmpg_uncc_edu::algorithms;
using namespace bmpg_uncc_edu::util::logger;

AtomPair::AtomPair(const string& name1,
		   const string& name2)
{
	if(name1.compare(name2) <= 0){
		atom1_ = name1;
		atom2_ = name2;
	} else {
		atom1_ = name2;
		atom2_ = name1;
	}
}


AtomPair::AtomPair(const char* name1,
		   const char* name2)
{
	string s1 = name1, s2 = name2;
	if(s1.compare(s2) <= 0){
		atom1_ = s1;
		atom2_ = s2;
	} else {
		atom1_ = s2;
		atom2_ = s1;
	}
}

bool AtomPair::operator==(const AtomPair& rhs) const
{
	return (atom1_ == rhs.atom1_) && (atom2_ == rhs.atom2_);
}


bool AtomPair::operator<(const AtomPair& rhs) const
{
	if(atom1_ < rhs.atom1_){
		return true;
	} else if ( (atom1_ == rhs.atom1_) && atom2_ < rhs.atom2_){
		return true;
	}
	return false;
	
}

ostream& operator<<(ostream& out, const AtomPair& rhs)
{
	out << rhs.atom1_ << " " << rhs.atom2_;
	return out;
}

AminoAcidLibrary* AminoAcidLibrary::instance_ = NULL;

AminoAcidLibrary::AminoAcidLibrary() : is_loaded_(false)
{

}

AminoAcidLibrary* AminoAcidLibrary::instance()
{
	if(instance_ == NULL){
		instance_ = new AminoAcidLibrary();
	}
	return instance_;
}

void AminoAcidLibrary::add(AminoAcid* aa)
{
	string name = aa->short_name();
	AminoAcid *ptr = amino_acids_[name];					//check if an entry is already there
	if(ptr != NULL)
		delete ptr;		
	amino_acids_[name] = aa;
}
	

/**
	List of all amino acid names in the library.
	@return vector of amino acid names
*/
std::vector<string> AminoAcidLibrary::residue_names() const
{
	std::vector<string> names;
	AminoAcidLibrary::amino_acids_map_t::const_iterator it;
	for(it = amino_acids_.begin(); it != amino_acids_.end(); it++){
		names.push_back(it->first);
	}
	return names;
}

/**
	List of all amino acids in the library.
	@return vector of amino acid pointers
*/
std::vector<const AminoAcid*> AminoAcidLibrary::amino_acids() const
{
	std::vector<const AminoAcid*> aas;
	AminoAcidLibrary::amino_acids_map_t::const_iterator it;
	for(it = amino_acids_.begin(); it != amino_acids_.end(); it++){
		aas.push_back(it->second);
	}
	return aas;
}


/**
	Find an amino acid by name; returns NULL if not found.
	@param amino acid name
	@return AminoAcid pointer
*/
const AminoAcid* AminoAcidLibrary::find_by_sname(const string & short_name) const
{
	AminoAcidLibrary::amino_acids_map_t::const_iterator it = amino_acids_.find(short_name);
	if(it == amino_acids_.end())
		return NULL;
	return it->second;
}

/**
	Find an amino acid by its long name.
	@param amino acid long name
	@return AminoAcid pointer
*/
const AminoAcid* AminoAcidLibrary::find_by_lname(const string & long_name) const
{
	const AminoAcid* aa;
	AminoAcidLibrary::amino_acids_map_t::const_iterator it;
	for(it = amino_acids_.begin(); it != amino_acids_.end(); it++){
		aa = it->second;
		if(long_name == aa->long_name()){
			return aa;
		}
	}
	return NULL;
}

/**
	Load the library from an input stream.
	@param input stream
	@return true if successful; false if not
*/
bool AminoAcidLibrary::load(const string& fname)
{
	if(is_loaded_)
		clear();
	AALibLoader loader(fname,*this);
	loader.load();
	AminoAcid* aa = NULL;
	AminoAcidLibrary::amino_acids_map_t::iterator it;
	for(it = amino_acids_.begin(); it != amino_acids_.end(); it++){
		aa = it->second;
		aa->setup_chemical_distance_map();
	}
	is_loaded_ = true;
	return true;
}

/**
	Ensures that the library is properly loaded. Throws an exception if not loaded.
*/
void AminoAcidLibrary::ensure_loaded() const
{
	if (!is_loaded_) 
		throw bmpg_uncc_edu::util::Exception("AminoAcidLibrary::ensure_loaded: amino acid library not loaded yet.", __FILE__, __LINE__);
}

/**
	Free the memories allocated by the libary class.
*/
void AminoAcidLibrary::clear()
{
	AminoAcidLibrary::amino_acids_map_t::iterator it;
	for(it = amino_acids_.begin(); it != amino_acids_.end(); it++){
		AminoAcid *aa = it->second;
		if(aa != NULL)delete aa;
		aa = NULL;
	}

	amino_acids_.clear();
	is_loaded_ = false;
}
/**
	Destructor.
*/
AminoAcidLibrary::~AminoAcidLibrary()
{
	clear();
}

/**
	Overloaded output stream operator <<.
*/
ostream& operator<<(ostream & out, const AminoAcidLibrary & rhs)
{
	AminoAcidLibrary::amino_acids_map_t::const_iterator it;
	for(it = rhs.amino_acids_.begin(); it != rhs.amino_acids_.end(); it++){
		out << (*(*it).second) << endl;
	}
	return out;
}


AminoAcid::AminoAcid()
{
	
}

AminoAcid::AminoAcid(const string & s) : short_name_(s)
{
	
}

void AminoAcid::add_connection(const string& atom1,
			       const string& atom2)
{
	graph_.add_edge(atom1, atom2,1);					//weight "1" is the chemical distance
}

void AminoAcid::add_hybridization(const string& atom,
				  const string& hybridization)
{
	hybridization_map_[atom].insert(hybridization);
}
	
/**
	Get the hybridization states of an atom.
	@param atom atom
	@return hybridization states
*/
AminoAcid::hybridization_states_t AminoAcid::hybridization_states(const string& atom) const
{
	name_to_hybridization_map_t::const_iterator it;
	it = hybridization_map_.find(atom);
	if(it == hybridization_map_.end()){
		throw Exception("Atom " + atom + " is not protonizable", __FILE__, __LINE__);
	}
	return it->second;
}

/**
	Get all the atoms that can be (de-)protonated.
	@return name of atoms that can be (de-)protonated
*/
AminoAcid::atom_names_t AminoAcid::protonation_atoms() const
{
	atom_names_t result;
	
	name_to_hybridization_map_t::const_iterator it;
	for(it = hybridization_map_.begin(); it != hybridization_map_.end(); it++){
		result.push_back(it->first);	
	}
	return result;
}
	
	
	
/**
	Set up the chemical distance map for the animo acid.
	@param graph the amino acid graph
*/
void AminoAcid::setup_chemical_distance_map()
{
	typedef graph_t::vertex_descriptor_t vertex_descriptor_t;
	typedef graph_t::edge_weight_t edge_weight_t;
	typedef graph_t::vertex_weight_map_t distance_map_t;


	distance_map_t dist;
	distance_map_t::iterator itd;

	atom_names_t atoms = atom_names();
	atom_names_t::const_iterator it;
	for(it = atoms.begin(); it != atoms.end(); it++){
		dist = graph_.Dijkstra(*it);
		for(itd = dist.begin(); itd != dist.end(); itd++){
			string name = itd->first;
			
			AtomPair key(*it,name);
			
			chemical_distance_map_[key] = itd->second;
		}
	}
}

void AminoAcid::print_chemical_distance_map(ostream& out) const
{
	chemical_distance_map_t::const_iterator it;
	for(it = chemical_distance_map_.begin(); it != chemical_distance_map_.end(); it++){
		out << it->first << " : " << it->second << endl;
	}
}

/**
	Get the closest distance of an atom to a certain element type.
	@param name atom name
	@param element element
	@return atom name; empty string if not found
*/
int AminoAcid::closest_distance(const string& name,
				const Element* element) const
{
	atom_names_t names = atom_names(element);
	int min = 10000;
	atom_names_t::const_iterator it;
	for(it = names.begin(); it != names.end(); it++){
		int d = chemical_distance(name, *it);
		if(d < min){
			min = d;
		}
	}
	return min;
}

/**
	Checks if an atom has a neighbor of a certain element type within a certain range.
	@param name atom name
	@param element element
	@param degree distance range
	@return atom name; empty string if not found
*/
string AminoAcid::closest_neighbor(const string& name,
				   const Element* element,
				   int degree) const
{
	string result = "";
	atom_names_t names = atom_names(element);
	int min = 1000000;
	
	atom_names_t::const_iterator it;
	for(it = names.begin(); it != names.end(); it++){
		int d = chemical_distance(name, *it);
		if(d < min){
			min = d;
			result = *it;
		}
	}
	if(min <= degree)
		return result;
	
	return "";
}

/**
	Get the chemical distance between two atoms.
	@param name1 atom 1
	@param name2 atom 2
	@return chemical distance
*/
int AminoAcid::chemical_distance(const string& name1,
				 const string& name2) const
{
	chemical_distance_map_t::const_iterator it;

	AtomPair key(name1,name2);
	it = chemical_distance_map_.find(key);
	
	if(it == chemical_distance_map_.end()){
		const char* alt1 = NULL, *alt2 = NULL;
		if(!has_atom(name1)){
			alt1 = PDBHelper::alternate_hydrogen_name(name1.c_str());
		} else {
			alt1 = name1.c_str();
		}
		
		if(!has_atom(name2)){
			alt2 = PDBHelper::alternate_hydrogen_name(name2.c_str());
		} else {
			alt2 = name2.c_str();	
		}
		
		if(alt1 != NULL && alt2 != NULL){
			key = AtomPair(alt1,alt2);
			it = chemical_distance_map_.find(key);
		}
		
		if(it == chemical_distance_map_.end()){
			string msg = "AminoAcid::chemical_distance: chemical distance between " + name1 + " and " + name2 + " in " + short_name_ + " is not found";
			Logger* logger = LoggerFactory::default_logger();
			logger->warn(msg);
			return 0;						//do not calculate force for unrecognized atoms
		}
	}
	return it->second;
}	

/**
	Get list of atom names in the amino acid.
	@return list of atom names
*/
AminoAcid::atom_names_t AminoAcid::atom_names() const
{
	atom_names_t result;
	graph_t::vertex_list_t v = graph_.vertices();
	result.insert(result.begin(),v.begin(), v.end());
	return result;
}


/**
	Get list of atom names in the amino acid.
	@return list of atom names
*/
AminoAcid::atom_names_t AminoAcid::atom_names(const Element* e) const
{
	atom_names_t names = atom_names();
	atom_names_t result;
	atom_names_t::const_iterator it;
	for(it = names.begin(); it != names.end(); it++){
		if( e == PDBHelper::atom_to_element(*it))
			result.push_back(*it);
	}
	return result;
}


/**
	Calculates the number of covalent bonds connecting two Carbon atoms.
	@param name1 atom name 1
	@param name2 atom name 2
	@return number of covalent bonds connecting the two atoms; 0 if not connected.
	//FIXME - the current implementation is not correct
*/
size_t AminoAcid::bond_number(const char* name1,
			      const char* name2) const
{
	const Element* e1 = PDBHelper::atom_to_element(name1);
	const Element* e2 = PDBHelper::atom_to_element(name2);
	PeriodicTable* pt = PeriodicTable::instance();
	const Element* C = pt->element("C");
	if(e1 != C || e2 != C)
		return 0;
	
	size_t n1 = number_of_neighbors(name1);
	size_t n2 = number_of_neighbors(name2);
	
	size_t number = (8 - (n1 + n2))/2 + 1;
	return number;
}

AminoAcid::neighbor_list_t AminoAcid::nearest_neighbors(const string& name) const
{
	return graph_.neighbors(name);
}
	
AminoAcid::neighbor_list_t AminoAcid::nearest_neighbors(const char* name) const
{
	string s = name;
	return graph_.neighbors(s);
}

size_t AminoAcid::number_of_neighbors(const string& name) const
{
	return graph_.number_of_neighbors(name);
}

size_t AminoAcid::number_of_neighbors(const char* name) const
{
	string s = name;
	return number_of_neighbors(s);
}


/**
	Test if two atoms are neighbors based on the record in the library.
	@param Atom name 1
	@param Atom name 2
	@return true if neighbors; false if not
*/
bool AminoAcid::is_neighbor(const string & name1,
			    const string & name2) const
{
	string s1 = name1;
	string s2 = name2;
	return chemical_distance(s1,s2) == 1;
}

/**
	Test if an atom is neighboring with an element.
	@param name atom anme
	@param element element
	@return true if neighbors; false if not
*/
bool AminoAcid::is_neighbor(const string& name,
			    const Element* element) const
{
	atom_names_t names = atom_names(element);
	atom_names_t::const_iterator it;
	for(it = names.begin(); it != names.end(); it++){
		if(chemical_distance(name,*it) == 1)
			return true;
	}
	return false;
}

/**
	Overloaded output stream operator <<.
*/
ostream& operator<<(ostream & out, const AminoAcid & rhs)
{
	out << "Amino acid: " <<  setw(3) << rhs.short_name() << endl;
	
	out << rhs.graph_ << endl;
	
	return out;
}
	

/**
	Find an atom by name in the amino acid.
	@param atom name string
	@return pointer to atom
*/
bool AminoAcid::has_atom(const string & atom_name) const
{
	return graph_.has(atom_name);
}

/**
	Destructor.
*/
AminoAcid::~AminoAcid()
{
	formula_.clear();
}


}	//namespace bmpg_uncc_edu::chemistry::library	
}	//namespace bmpg_uncc_edu::chemistry
}	//namespace bmpg_uncc_edu

#endif

