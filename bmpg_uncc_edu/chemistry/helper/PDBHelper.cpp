/**
   Copyright (c) 2008 by Hui Wang
   @author Hui Wang
   @date Jan 20, 2009
*/


#ifndef bmpg_uncc_edu_chemistry_helper_PDBHelper_cpp
#define bmpg_uncc_edu_chemistry_helper_PDBHelper_cpp

#include <cstring>
#include <iostream>
#include <sstream>
#include <bmpg_uncc_edu/chemistry/helper/PDBHelper.hpp>
#include <bmpg_uncc_edu/chemistry/helper/TorsionAngle.hpp>
#include <bmpg_uncc_edu/chemistry/Element.hpp>
#include <bmpg_uncc_edu/chemistry/PDBAtom.hpp>
#include <bmpg_uncc_edu/chemistry/PDBProteinResidue.hpp>
#include <bmpg_uncc_edu/chemistry/PDBProtein.hpp>
#include <bmpg_uncc_edu/chemistry/library/PeriodicTable.hpp>
#include <bmpg_uncc_edu/util/logger/LoggerFactory.hpp>
#include <bmpg_uncc_edu/math/MathConstants.hpp>
#include <bmpg_uncc_edu/util/Exception.hpp>

namespace bmpg_uncc_edu {
namespace chemistry {
namespace helper {
using namespace bmpg_uncc_edu::math;
using namespace bmpg_uncc_edu::chemistry::library;
using namespace bmpg_uncc_edu::util::logger;
	
const Element* PDBHelper::H_ = NULL;
	
const char * PDBHelper::alternate_hydrogen_name(const char * name)
{
	size_t i;
	const size_t num_pairs = 39;
	const char * alt_h_names[num_pairs][2] = {
		{" HA1","1HA "}, {" HA2","2HA "}, {" HA3","3HA "}, 
		{" HB1","1HB "}, {" HB2","2HB "}, {" HB3","3HB "},
		{" HG1","1HG "}, {" HG2","2HG "}, {" HG3","3HG "},
		{" HD1","1HD "}, {" HD2","2HD "}, {" HD3","3HD "},
		{" HE1","1HE "}, {" HE2","2HE "}, {" HE3","3HE "},
		{" HZ1","1HZ "}, {" HZ2","2HZ "}, {" HZ3","3HZ "},
		{" HH1","1HH "}, {" HH2","2HH "}, {" HH3","3HH "},
		{"HG11","1HG1"}, {"HG12","2HG1"}, {"HG13","3HG1"},
		{"HG21","1HG2"}, {"HG22","2HG2"}, {"HG23","3HG2"},
		{"HD11","1HD1"}, {"HD12","2HD1"}, {"HD13","3HD1"},
		{"HD21","1HD2"}, {"HD22","2HD2"}, {"HD23","3HD2"},
		{"HE21","1HE2"}, {"HE22","2HE2"},
		{"HH11","1HH1"}, {"HH12","2HH1"},
		{"HH21","1HH2"}, {"HH22","2HH2"}
	};
	
	for (i = 0; i < num_pairs; ++i) {
		if (strncmp(alt_h_names[i][0], name, 4) == 0){
			return alt_h_names[i][1];
		}
		if (strncmp(alt_h_names[i][1], name, 4) == 0){
			return alt_h_names[i][0];
		}
	}
	
	return NULL;
}

const char * PDBHelper::alternate_hydrogen_name(const string& name)
{
	if(name.size() < 4){
		return NULL;	
	}
	return alternate_hydrogen_name(name.c_str());
}

string PDBHelper::alternate_his_name(Handle<PDBProteinResidue> res)
{
	if(res.empty())
		throw bmpg_uncc_edu::util::Exception("PDBHelper::alternate_his_name: NULL residue pointer", __FILE__, __LINE__);
	
	if(res->aa->short_name() != "HIS")
		cerr << "PDBHelper::alternate_his_name: not a HIS residue\n";
	
	bool has_hd1 = false, has_he2 = false;
	
	PDBProteinResidue::atoms_t::const_iterator it;
	for(it = res->atoms().begin(); it != res->atoms().end(); it++){
		string name = (*it)->atom_name;
		if(name == " HE2"){
			has_he2 = true;
			continue;
		} else if(name == " HD1"){
			has_hd1 = true;
			continue;
		}
	}
	
	if(has_hd1 && has_he2)
		return "HIS";
	else if(has_hd1)
		return "HID";
	else if(has_he2)
		return "HIE";
	else {
		return "";							//missing both hD1 and HE2
	}
}

/**
	Use the first two letters of the four-letter PDBAtom name to search for the corresponding element.
	First check if the name starts with a number or space, if so, use the second letter only;
	Second, check if the second letter is a number, if so, use the first letter only; e.g., this could happen for H0A;
	Third, check if the four-letter name contains number AND 'H', if so, it is H;
	Fourth, if the element is not found, check if the name contains 'H', if so, it is likely to be a hydrogen; e.g., 'HA'
	If all of the above tests fail, then the element name is not recognizable, and an exception is thrown.
	@param aname atom name
	@return element
*/
const Element* PDBHelper::atom_to_element(const string& sname)
{
	library::PeriodicTable *pt = library::PeriodicTable::instance();	
	string s = sname.substr(0,2);						//get the first two characters
	
	if(s.find_first_of(" 0123456789") == 0){				//first character is a number of space, use the second character only
		s = s.substr(1,1);
	} else if(s.find_first_of("0123456789") == 1){
		s = s.substr(0,1);						//H0 is an H
	}
	
	if(sname.find_first_of("0123456789") != string::npos &&	
		s.find_first_of("Hh") != string::npos){				//contains number and H
		return pt->element("H");
	}
		
	const Element* e = pt->element(s);
	if(e != NULL){
		return e;
	} else if(s.find_first_of("hH") != string::npos){			//found H
		return pt->element("H");
	} else {								//no H
		string msg = "Unrecognized element " + s + "; returning NULL";
		Logger *logger = LoggerFactory::default_logger();
		logger->critical(msg);
	}
	return NULL;
}

const Element* PDBHelper::atom_to_element(const char* aname)
{
	string sname(aname);
	return atom_to_element(sname);	
}

const Element* PDBHelper::atom_to_element(const PDBAtom* atom)
{
	if(atom == NULL)
		throw Exception("NULL pointer for atom.", __FILE__, __LINE__);
	
	string sname(atom->atom_name);
	return atom_to_element(sname);	
}

bool PDBHelper::is_c_alpha(const PDBAtom* atom)
{
	string name = atom->atom_name;
	return (name.compare(" CA ") == 0);
}

bool PDBHelper::is_heavy_atom(const PDBAtom* atom)
{
	if(H_ == NULL){
		library::PeriodicTable* pt = library::PeriodicTable::instance();
		H_ = pt->element("H");
	}
	return (atom->element() != H_);
}

bool PDBHelper::is_heavy_atom(const string& atom)
{
	if(H_ == NULL){
		library::PeriodicTable* pt = library::PeriodicTable::instance();
		H_ = pt->element("H");
	}
	return (atom_to_element(atom) != H_);
}

bool PDBHelper::is_hydrogen_atom(const PDBAtom* atom)
{
	return !is_heavy_atom(atom);	
}
/**
	Check if an atom is on the backbone.
	@param atom atom
	@return 0 if not; 1 if 'N', 2 if "CA", 3 if 'C'
*/
int PDBHelper::is_back_bone(const PDBAtom* atom)
{
	if(atom == NULL)
		throw bmpg_uncc_edu::util::Exception("PDBHelper::is_back_bone: atom is NULL.", __FILE__, __LINE__);
	
	const AminoAcid* aa = atom->aa;
	
	if( strcmp(atom->atom_name," CA ") == 0
		&& aa->chemical_distance(atom->atom_name, " N  ") == 1
		&& aa->chemical_distance(atom->atom_name, " C  ") == 1
	)
	return 2;
	
	if( strcmp(atom->atom_name," N  ") == 0
		&& aa->chemical_distance(atom->atom_name, " CA ") == 1
		&& aa->chemical_distance(atom->atom_name, " C  ") <= 2
	)
	return 1;
	
	
	if( strcmp(atom->atom_name," C  ") == 0
		&& aa->chemical_distance(atom->atom_name, " CA ") == 1
		&& aa->chemical_distance(atom->atom_name, " N  ") == 2
	)
	return 3;
	
	return 0;
}

/**
	Checks if an atom is on the backbone, including carbonyl oxygen.
	@param atom atom
	@return 0 if not on the main chain; positive integer if is on the main chain.
*/
int PDBHelper::is_main_chain(const PDBAtom* atom){

	const AminoAcid* aa = atom->aa;
	
	int backbone = PDBHelper::is_back_bone(atom);
  
	if( backbone )
		return( backbone );
	else if(strcmp(atom->atom_name," O  ") == 0 && 
	   aa->chemical_distance(atom->atom_name, " C  ") == 1 &&
	   aa->chemical_distance(atom->atom_name, " CA ") <= 2)
	return(4);
	else if(strcmp(atom->atom_name," H  ") == 0 &&
	   aa->chemical_distance(atom->atom_name, " N  ") == 1 &&
	   aa->chemical_distance(atom->atom_name, " CA ") <= 2)
	return(5);
  return(0);
}

bool PDBHelper::is_C_terminal(Handle<PDBProteinResidue> res)
{
	return (res->next_res).empty();
}


bool PDBHelper::is_N_terminal(Handle<PDBProteinResidue> res)
{
	return (res->prev_res).empty();
}

bool PDBHelper::is_terminal(Handle<PDBProteinResidue> res)
{
	return (res->prev_res).empty() || (res->next_res).empty();
}

double PDBHelper::distance(const PDBAtom* atom1,
			   const PDBAtom* atom2)
{
	double dx = atom1->x - atom2->x;
	double dy = atom1->y - atom2->y;
	double dz = atom1->z - atom2->z;
	return sqrt(dx*dx + dy*dy + dz*dz);
}


string PDBHelper::residue_info(Handle<PDBProteinResidue> res)
{
	return residue_info(*res);
}


string PDBHelper::residue_info(const PDBProteinResidue& res)
{
	stringstream ss;
	string filename = res.file->filename;
	if(res.insertion_code == ' '){
		ss << "[" << filename << "---" << res.number << "--" << res.aa->short_name() << "]";
	} else {
		ss << "[" << filename << "-" << res.insertion_code << "-" << res.number << "]";	
	}
	return ss.str();	
}
	
bool PDBHelper::get_atoms(std::vector<const PDBAtom*>& atoms,
			  Handle<PDBProteinResidue> res,
			  const std::vector<std::string> atom_names,
			  const char * res_ptrs,
			  bool try_alt_hydrogen_names)
{
	size_t n = atoms.size();
	
	if (strlen(res_ptrs) != n){
		stringstream ss;
		ss << "PDBHelper::get_atoms(): res_ptrs must have exactly " << n << " characters specified";
		throw bmpg_uncc_edu::util::Exception(ss.str(), __FILE__, __LINE__);
	}
	
	Handle<PDBProteinResidue> res_ptr;
	size_t i = 0;
	for (i = 0; i < n; ++i) {
		if (res_ptrs[i] == 'P') {
			if ((res_ptr = res->prev_res).empty())
				return false;
		}
		else if (res_ptrs[i] == 'C')
			res_ptr = res;
		else if (res_ptrs[i] == 'N') {
			if ((res_ptr = res->next_res).empty())
				return false;
		}else if(res_ptrs[i] == 'S') {
			if ((res_ptr = res->sulf_res).empty())
				return false;
		}
		else
			throw bmpg_uncc_edu::util::Exception("PDBHelper::get_atoms(): invalid res_ptr: must be 'N', 'C', 'P', or 'S'.", __FILE__, __LINE__);
	
		
		if ((atoms[i] = res_ptr->atom(atom_names[i].c_str(),try_alt_hydrogen_names)) == NULL)
			return false;
		
	}
	return true;
}

/**
	Calculate the torsion angle formed by four atoms.
	@param atoms list of atoms
	@param tolerance
	@return angle in range [-pi,pi]	
*/
double PDBHelper::torsion_angle(const std::vector<const PDBAtom*>& atoms,
				double tol)
{
	if (atoms[0] == NULL || atoms[1] == NULL
		|| atoms[2] == NULL || atoms[3] == NULL) 
	{
		throw bmpg_uncc_edu::util::Exception("PDBHelper::torsion_angle(): NULL PDBAtom specified.", __FILE__, __LINE__);
	}
	
	typedef bmpg_uncc_edu::math::R3Vector<double> vec_t;
	vec_t v0(atoms[0]->x, atoms[0]->y, atoms[0]->z);
	vec_t v1(atoms[1]->x, atoms[1]->y, atoms[1]->z);
	vec_t v2(atoms[2]->x, atoms[2]->y, atoms[2]->z);
	vec_t v3(atoms[3]->x, atoms[3]->y, atoms[3]->z);
	
	return torsion_angle(v0,v1,v2,v3,tol);
}


double PDBHelper::torsion_angle(double atom1x,
				double atom1y,
				double atom1z,
				double atom2x,
				double atom2y,
				double atom2z,
				double atom3x,
				double atom3y,
				double atom3z,
				double atom4x,
				double atom4y,
				double atom4z,
				double tol)
{
	typedef bmpg_uncc_edu::math::R3Vector<double> vec_t;
	
	// Setup position vectors for each of the four atom
	vec_t v1(atom1x, atom1y, atom1z);
	vec_t v2(atom2x, atom2y, atom2z);
	vec_t v3(atom3x, atom3y, atom3z);
	vec_t v4(atom4x, atom4y, atom4z);
	// Call the vector version of this function
	return torsion_angle(v1, v2, v3, v4, tol);
}

double PDBHelper::torsion_angle(const bmpg_uncc_edu::math::R3Vector<double> & v1,
				const bmpg_uncc_edu::math::R3Vector<double> & v2,
				const bmpg_uncc_edu::math::R3Vector<double> & v3,
				const bmpg_uncc_edu::math::R3Vector<double> & v4,
				double tol)
{
	typedef bmpg_uncc_edu::math::R3Vector<double> vec_t;
	
	double acos_argument, angle;
	vec_t v12(v2-v1);
	vec_t v23(v3-v2);
	vec_t v34(v4-v3);
	
	vec_t n12(v12.cross(v23).normalize());
	vec_t n34(v23.cross(v34).normalize());
	
	// Why do I (mjf) not directly call n12.angle_to(n34)? 
	// Answer: Because of poor precision in the atomic coordinates,
	// the argument to std::acos() might be slightly larger than 1 (e.g. 1.0001)
	// or slightly smaller than -1 (e.g. -1.0001).  When this happens, the
	// call to std::acos() returns NAN.  Therefore, I handle this calculation
	// explicitly, allowing a user-specified tolerance.  If the argument to
	// the std::acos() function is more than tol larger (smaller) than 1
	// (or -1), then NAN is returned anyway.  The tolerance provides a region
	// of "expected" roundoff error on account of poor precision in atomic
	// coordinates.
	acos_argument = (n12.inner_product(n34))/(n12.norm()*n34.norm());
	if (acos_argument < -1 && acos_argument > -(1+tol)) 
		angle = std::acos(-1.0);
	else if (acos_argument > 1 && acos_argument < (1+tol))
		angle = std::acos(1.0);
	else
		angle = std::acos(acos_argument);

	// Possibly update the sign of the result
	vec_t n12_n34 = n12.cross(n34);
	if (n12_n34.inner_product(v23) < 0.0)
		angle = -angle;
	return angle;
}


/**
	Calculate the out-of-plane angle formed by four atoms.
	@param atoms list of atoms
	@param tolerance
	@return angle in range [-pi,pi]	
*/
double PDBHelper::oop_angle(const std::vector<const PDBAtom*>& atoms,
			    double tol)
{
	if (atoms[0] == NULL || atoms[1] == NULL
		|| atoms[2] == NULL || atoms[3] == NULL) 
	{
		throw bmpg_uncc_edu::util::Exception("PDBHelper::oop_angle(): NULL PDBAtom specified.", __FILE__, __LINE__);
	}
	
	typedef bmpg_uncc_edu::math::R3Vector<double> vec_t;
	vec_t v0(atoms[0]->x, atoms[0]->y, atoms[0]->z);
	vec_t v1(atoms[1]->x, atoms[1]->y, atoms[1]->z);
	vec_t v2(atoms[2]->x, atoms[2]->y, atoms[2]->z);
	vec_t v3(atoms[3]->x, atoms[3]->y, atoms[3]->z);
	
	return oop_angle(v0,v1,v2,v3,tol);
}

double PDBHelper::oop_angle(const bmpg_uncc_edu::math::R3Vector<double> & v1,
			    const bmpg_uncc_edu::math::R3Vector<double> & v2,
			    const bmpg_uncc_edu::math::R3Vector<double> & v3,
			    const bmpg_uncc_edu::math::R3Vector<double> & v4,
			    double tol)
{
	typedef bmpg_uncc_edu::math::R3Vector<double> vec_t;
	
	double acos_argument, angle;
	vec_t v13(v3-v1);
	vec_t v23(v3-v2);
	vec_t v34(v4-v3);
	
	vec_t n12(v13.cross(v23).normalize());
	vec_t n34(v34.normalize());
	
	acos_argument = (n12.inner_product(n34))/(n12.norm()*n34.norm());
	if (acos_argument < -1 && acos_argument > -(1+tol)) 
		angle = std::acos(-1.0);
	else if (acos_argument > 1 && acos_argument < (1+tol))
		angle = std::acos(1.0);
	else
		angle = std::acos(acos_argument);

	angle = abs(angle - MathConstants::pi/2);

	return angle;
}


/**
	Calculate the distance between two atoms.
	@param atoms list of atoms
	@return distance
*/	
double PDBHelper::bond_length(const std::vector<const PDBAtom*>& atoms)
{
	if (atoms[0] == NULL || atoms[1] == NULL) 
	{
		throw bmpg_uncc_edu::util::Exception("PDBHelper::bond_length(): NULL PDBAtom specified.", __FILE__, __LINE__);
	}
	
	typedef bmpg_uncc_edu::math::R3Vector<double> vec_t;
	vec_t v0(atoms[0]->x, atoms[0]->y, atoms[0]->z);
	vec_t v1(atoms[1]->x, atoms[1]->y, atoms[1]->z);
	
	vec_t v10(v0 - v1);
	double length = v10.norm();
	return length;
}

/**
	Calculate the angle formed by three atoms.
	@param atoms list of atoms
	@param tol tolerance
	@return angle in range [0,pi]
*/	
double PDBHelper::bending_angle(const std::vector<const PDBAtom*>& atoms)
{
	if (atoms[0] == NULL || atoms[1] == NULL
		|| atoms[2] == NULL) 
	{
		throw bmpg_uncc_edu::util::Exception("PDBHelper::bending_angle(): NULL PDBAtom specified.", __FILE__, __LINE__);
	}
	
	typedef bmpg_uncc_edu::math::R3Vector<double> vec_t;
	vec_t v0(atoms[0]->x, atoms[0]->y, atoms[0]->z);
	vec_t v1(atoms[1]->x, atoms[1]->y, atoms[1]->z);
	vec_t v2(atoms[2]->x, atoms[2]->y, atoms[2]->z);
	
	vec_t v10(v0 - v1);
	vec_t v12(v2 - v1);
	v10.normalize();
	v12.normalize();
	
	double dot = v10.inner_product(v12);
	double angle = std::acos(dot);
	return angle;
}


/**
	Calculate the angle formed by three atoms.
	@param atoms list of atoms
	@param tol tolerance
	@return angle in range [0,pi]
*/	
double PDBHelper::bending_angle(const bmpg_uncc_edu::math::R3Vector<double> & v0,
				const bmpg_uncc_edu::math::R3Vector<double> & v1,
				const bmpg_uncc_edu::math::R3Vector<double> & v2)
{	
	
	bmpg_uncc_edu::math::R3Vector<double> v10(v0 - v1);
	bmpg_uncc_edu::math::R3Vector<double> v12(v2 - v1);
	v10.normalize();
	v12.normalize();
	
	double dot = v10.inner_product(v12);
	double angle = std::acos(dot);
	return angle;
}

/**
	Calculate the angle formed by three atoms.
	@param atom1 first atom
	@param atom2 central atom
	@param atom3 thrid ataom
	@return angle in range [0,pi]
*/	
double PDBHelper::bending_angle(const PDBAtom* atom1,
				const PDBAtom* atom2,
				const PDBAtom* atom3)
{
	if (atom1 == NULL || atom2 == NULL
		|| atom3 == NULL) 
	{
		throw bmpg_uncc_edu::util::Exception("PDBHelper::bending_angle(): NULL PDBAtom specified.", __FILE__, __LINE__);
	}
	
	std::vector<const PDBAtom*> atoms;
	atoms.push_back(atom1);
	atoms.push_back(atom2);
	atoms.push_back(atom3);
	return bending_angle(atoms);
}

double PDBHelper::torsion_angle_phi(Handle<PDBProteinResidue> res)
{
	std::vector<const PDBAtom*> atoms;
	atoms.push_back(NULL);
	atoms.push_back(NULL);
	atoms.push_back(NULL);
	atoms.push_back(NULL);
	
	std::vector<std::string> atom_names;						//the names of the four atoms involved
	const char * res_ptrs = "PCCC";						//pointers to the residues that each atom belongs to
	atom_names.push_back(" C  ");
	atom_names.push_back(" N  ");
	atom_names.push_back(" CA ");
	atom_names.push_back(" C  ");
	
	bool no_missing_atom_flag = PDBHelper::get_atoms(atoms, res, atom_names, res_ptrs, false);
	
	if(!no_missing_atom_flag) return 0;					//FIXME - deal with N and C terminals
	return PDBHelper::torsion_angle(atoms,0.01);
}

double PDBHelper::torsion_angle_psi(Handle<PDBProteinResidue> res)
{
	std::vector<const PDBAtom*> atoms;
	atoms.push_back(NULL);
	atoms.push_back(NULL);
	atoms.push_back(NULL);
	atoms.push_back(NULL);
	
	std::vector<std::string> atom_names;						//the names of the four atoms involved
	const char * res_ptrs = "CCCN";						//pointers to the residues that each atom belongs to
	atom_names.push_back(" N  ");
	atom_names.push_back(" CA ");
	atom_names.push_back(" C  ");
	atom_names.push_back(" N  ");
	
	
	bool no_missing_atom_flag = PDBHelper::get_atoms(atoms, res, atom_names, res_ptrs, false);
	
	if(!no_missing_atom_flag) return 0;					//FIXME - deal with N and C terminals
	return PDBHelper::torsion_angle(atoms,0.01);
}


double PDBHelper::torsion_angle_omega(Handle<PDBProteinResidue> res)
{
	std::vector<const PDBAtom*> atoms;
	atoms.push_back(NULL);
	atoms.push_back(NULL);
	atoms.push_back(NULL);
	atoms.push_back(NULL);
	
	std::vector<std::string> atom_names;						//the names of the four atoms involved
	const char * res_ptrs = "CCNN";						//pointers to the residues that each atom belongs to
	atom_names.push_back(" CA ");
	atom_names.push_back(" C  ");
	atom_names.push_back(" N  ");
	atom_names.push_back(" CA ");
	
	bool no_missing_atom_flag = PDBHelper::get_atoms(atoms, res, atom_names, res_ptrs, false);
	
	if(!no_missing_atom_flag) return 0;					//FIXME - deal with N and C terminals
	return PDBHelper::torsion_angle(atoms,0.01);
}

void PDBHelper::print_xyz(const PDBProtein& protein, std::ostream& os)
{
	os << protein.num_atoms() << "\n";
	PDBProtein::atom_iterator_t ita;
	for (ita = protein.atom_begin(); ita != protein.atom_end(); ++ita) {
		const PDBAtom* a = *ita;
		os << a->x << "\t" << a->y << "\t" << a->z << "\t" << a->atom_name << "\n";
	}
}

} // namespace bmpg_uncc_edu::chemistry::helper
} // namespace bmpg_uncc_edu::chemistry
} // namespace bmpg_uncc_edu

#endif

