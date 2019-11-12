/**
  Copyright (c) 2008 by Mike Fairchild
  @author: Mike Fairchild, http://www.mikef.org, mjfairch@gmail.com
*/

#ifndef bmpg_uncc_edu_chemistry_TorsionAngle_hpp
#define bmpg_uncc_edu_chemistry_TorsionAngle_hpp

#include <string>
#include <map>
#include <vector>


namespace bmpg_uncc_edu {
	namespace chemistry {
		namespace library {
			class AminoAcid;
		}
	}
}

namespace bmpg_uncc_edu {
namespace chemistry {
namespace helper {
using namespace bmpg_uncc_edu::chemistry::library;

struct TorsionAngleSpecification {
	std::string pdb_atom_names[4];
	char res_ptr[5];
	
	TorsionAngleSpecification() { clear(); }
	~TorsionAngleSpecification() { clear(); }
	
	void clear() {
		pdb_atom_names[0] = pdb_atom_names[1] = "";
		pdb_atom_names[2] = pdb_atom_names[3] = "";
		res_ptr[0] = res_ptr[1] = res_ptr[2] = res_ptr[3] = 'C';
	}
};

// This is a Singleton class, like the AminoAcidDatabase.
// See: http://en.wikipedia.org/wiki/Singleton_pattern
//
class TorsionAngleDatabase
{	
public:
	// Error codes returned on load() method.
	static const int ERR_LINE_TOO_SHORT;		// Defined in .cpp file
	static const int ERR_INVALID_AMINO_ACID;	// Defined in .cpp file
	static const int ERR_INVALID_RES_PTR;		// Defined in .cpp file
	
	~TorsionAngleDatabase();
	static TorsionAngleDatabase * instance();
	
	inline bool is_loaded() const { return is_loaded_; }
	int load(std::istream & is);
	void clear();
	
	bool has_angle(const AminoAcid * aa,
		       const std::string & angle_name) const;
	
	bool has_angle(const std::string & res_sname,
		       const std::string & angle_name) const;
	
	bool get_angle_names(const AminoAcid * aa,
			     std::vector<std::string> & angle_names) const;
	
	bool get_angle_names(const std::string & res_sname,
			     std::vector<std::string> & angle_names) const;
	
	bool get_torsion_angle_spec(const AminoAcid * aa,
				    const std::string & angle_name,
				    const TorsionAngleSpecification * & ta_spec) const;
	
	bool get_torsion_angle_spec(const std::string & res_sname,
				    const std::string & angle_name,
				    const TorsionAngleSpecification * & ta_spec) const;

protected:
	inline void ensure_loaded() const;		// Throws exception if is_loaded_ == false
	inline void ensure_aa_loaded() const;	// Throws exception if amino acid datbase not loaded

private:
	// Hide constructors and assignment operator to ensure Singleton pattern
	TorsionAngleDatabase();
	TorsionAngleDatabase(const TorsionAngleDatabase & rhs);
	TorsionAngleDatabase & operator=(const TorsionAngleDatabase & rhs);
	
	static TorsionAngleDatabase * instance_;
	typedef std::map<const std::string, TorsionAngleSpecification *> ta_set_t;
	typedef std::map<const AminoAcid *, ta_set_t *> ta_map_t;
	
	ta_map_t ta_map_;
	bool is_loaded_;
};


} // namespace bmpg_uncc_edu::chemistry::helper
} // namespace bmpg_uncc_edu::chemistry
} // namespace bmpg_uncc_edu

#endif

