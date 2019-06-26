/**
   Copyright (c) 2008 by Hui Wang
   @author Hui Wang
   @date Jan 20, 2009
*/


#ifndef bmpg_uncc_edu_chemistry_library_AminoAcidLibrary_hpp
#define bmpg_uncc_edu_chemistry_library_AminoAcidLibrary_hpp

#include <set>
#include <vector>
#include <map>
#include <bmpg_uncc_edu/algorithms/Graph.hpp>

namespace bmpg_uncc_edu {
	namespace chemistry {
		class Element;
	}
}

namespace bmpg_uncc_edu {
namespace chemistry {
namespace library {
using namespace std;
using namespace bmpg_uncc_edu::algorithms;


/**
 * AtomPair defines two atoms linked chemically and is used to define their distance.
 */

class AtomPair
{
public:
	AtomPair(const string& name1,
		 	 const string& name2);
	
	AtomPair(const char* name1,
		 const char* name2);
	
	bool operator==(const AtomPair& rhs) const;
	bool operator<(const AtomPair& rhs) const;
	friend ostream& operator<<(ostream& out, const AtomPair& rhs);
	
private:
	string atom1_, atom2_;
};

/**
 * AminoAcid stores all the characteristics and properties of a single standard amino acid.
 */

class AminoAcid
{
public:
	typedef enum {acid, base, neutral} acid_base_t;
	
	typedef std::vector<string> atom_names_t;
	typedef std::set<string> hybridization_states_t;
	
	typedef std::map<AtomPair, int> chemical_distance_map_t;
	typedef std::map<string, hybridization_states_t> name_to_hybridization_map_t;
	
	typedef Graph<string,int> graph_t;
	typedef graph_t::vertex_list_t neighbor_list_t;
	
	char letter;								// Standard one character symbol
	bool polar;								// True if residue is polar, False otherwise
	acid_base_t acid_base;								// A=acid, B=base, or N=neutral
	bool branched;								// True if residue has a branch (ring or tree)
	
	
	AminoAcid();
	AminoAcid(const string & s);
	
	void add_connection(const string& atom1,
			    const string& atom2);
	
	void add_hybridization(const string& atom,
			       const string& hybridization);
	
	hybridization_states_t hybridization_states(const string& atom) const;
	atom_names_t protonation_atoms() const;
	
	
	void setup_chemical_distance_map();

	atom_names_t atom_names() const;
	atom_names_t atom_names(const Element* e) const;
	
	string closest_neighbor(const string& name,
				const Element* element,
				int degree) const;
	
	int closest_distance(const string& name,
			     const Element* element) const;
	
	size_t bond_number(const char* name1,
			   const char* name2) const;

	bool has_atom(const string & atom_name) const;
	
	bool is_neighbor(const string & name1,
			 const string & name2) const;
	
	bool is_neighbor(const string & name1,
			 const Element* element) const;
	
	bool is_acid() const{return acid_base == acid;}
	bool is_base() const{return acid_base == base;}
	bool is_neutral() const{return acid_base == neutral;}
	
	neighbor_list_t nearest_neighbors(const string& name) const;
	neighbor_list_t nearest_neighbors(const char* name) const;
	size_t number_of_neighbors(const string& name) const;
	size_t number_of_neighbors(const char* name) const;
	string& formula() {return formula_;}
	string& short_name() {return short_name_;}
	string& long_name() {return long_name_;}
	string long_name() const {return long_name_;}
	string short_name() const {return short_name_;}
	
	double naked_asa() const {return naked_asa_;}
	double& naked_asa() {return naked_asa_;}
	
	double flanking_asa() const {return flanking_asa_;}
	double& flanking_asa() {return flanking_asa_;}
	
	
	int chemical_distance(const string& name1,
			      const string& name2) const;			//calculate the chemical distance
	
	void print_chemical_distance_map(ostream& out) const;

	graph_t graph() const{return graph_;}
	
	friend ostream& operator<<(ostream & out, const AminoAcid & rhs);
	virtual ~AminoAcid();

private:
	string short_name_, long_name_;
	chemical_distance_map_t chemical_distance_map_;
	string formula_;
	graph_t graph_;
	name_to_hybridization_map_t hybridization_map_;
	double naked_asa_, flanking_asa_;
};



/**
 * AminoAcidLibrary stores the standard amino acids. Each amino acid stores the atom labels, and the covalent bonds
 * among the atoms in the same amino acid. PeriodicTable must be loaded before loading this class.
*/
class AminoAcidLibrary 
{
public:
	
	typedef std::map<string,AminoAcid*> amino_acids_map_t;
	typedef amino_acids_map_t::iterator iterator;
	typedef std::set<string> amino_acid_names_t;
	
	amino_acids_map_t amino_acids_;				//keyed by short name; stored in a map for quick searching
	
	
	iterator begin() {return amino_acids_.begin();}
	iterator end() {return amino_acids_.end();}

	static AminoAcidLibrary* instance();

	bool load(const string& fname);
	void add(AminoAcid* aa);
	
	void ensure_loaded() const;
	bool is_loaded() const {return is_loaded_;}
	size_t size() const {return amino_acids_.size();}
	bool has_amino_acid(const string & name) const {return amino_acids_.find(name) != amino_acids_.end();}
	
	std::vector<string> residue_names() const;
	std::vector<const AminoAcid*> amino_acids() const;

	const AminoAcid* find_by_lname(const string & name) const;
	const AminoAcid* find_by_sname(const string & short_name) const;
	bool has(const string& short_name) const{return find_by_sname(short_name) != NULL;}
	
	void add_amino_acid_to_fed(const string& aa_name) {amino_acid_names_.insert(aa_name);}
	
	amino_acid_names_t amino_acids_to_fed() const{return amino_acid_names_;}
	
	void clear();
	friend ostream& operator<<(ostream & out, const AminoAcidLibrary & rhs);
	virtual ~AminoAcidLibrary();
private:
	bool is_loaded_;
	static AminoAcidLibrary* instance_;
	amino_acid_names_t amino_acid_names_;
private:
	AminoAcidLibrary();						//hide default constructor
	AminoAcidLibrary(const AminoAcidLibrary &rhs);			//hide copy constructor
	AminoAcidLibrary& operator=(const AminoAcidLibrary &rhs);	//hide assignment operator
};


}	//namespace bmpg_uncc_edu::chemistry::library
}	//namespace bmpg_uncc_edu::chemistry
}	//namespace bmpg_uncc_edu

#endif

