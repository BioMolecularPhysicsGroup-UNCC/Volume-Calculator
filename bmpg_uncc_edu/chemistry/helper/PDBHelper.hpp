/**
   Copyright (c) 2008 by Hui Wang
   @author Hui Wang
   @date Jan 20, 2009
*/


#ifndef bmpg_uncc_edu_chemistry_PDBHelper_hpp
#define bmpg_uncc_edu_chemistry_PDBHelper_hpp

#include <vector>
#include <bmpg_uncc_edu/util/Handle.hpp>
#include <bmpg_uncc_edu/math/R3Vector.hpp>


namespace bmpg_uncc_edu {
	namespace chemistry {
		class PDBAtom;
		class Element;
		class PDBProteinResidue;
		class PDBProtein;
	}
}

namespace bmpg_uncc_edu {
namespace chemistry {
namespace helper {
using namespace std;
using namespace bmpg_uncc_edu::util;
using namespace bmpg_uncc_edu::math;

/**
 * PDBHelper is a collection of static functions useful for PDBProtein objects.  Examples include bond information, standard  
 * angles of various types, as well as identification of the presence of the existence of things such as C-alpha or backbone 
 * elements in a PDBAtom.
 */

class PDBHelper
{
public:
	static const char * alternate_hydrogen_name(const char * name);
	static const char * alternate_hydrogen_name(const string& name);
	
	static bool get_atoms(std::vector<const PDBAtom*>& atoms,
			      Handle<PDBProteinResidue> res,
			      const std::vector<std::string> atom_names,
			      const char * res_ptrs,
			      bool try_alt_hydrogen_names = true);
	
	static string alternate_his_name(Handle<PDBProteinResidue> res);
	static const Element* atom_to_element(const PDBAtom* atom);
	static const Element* atom_to_element(const char* aname);
	static const Element* atom_to_element(const string& aname);

	static double distance(const PDBAtom* aotm1,
			       const PDBAtom* atom2);
	
	static bool is_c_alpha(const PDBAtom* atom);
	static int is_back_bone(const PDBAtom* atom);
	static int is_main_chain(const PDBAtom* atom);
	static bool is_N_terminal(Handle<PDBProteinResidue> res);
	static bool is_C_terminal(Handle<PDBProteinResidue> res);
	static bool is_terminal(Handle<PDBProteinResidue> res);
	static bool is_heavy_atom(const PDBAtom* atom);
	static bool is_hydrogen_atom(const PDBAtom* atom);	
	static bool is_heavy_atom(const string& atom);

	
	static string residue_info(Handle<PDBProteinResidue> res);
	static string residue_info(const PDBProteinResidue& res);

	static void print_xyz(const PDBProtein& protein,
			      std::ostream& os);
	
	static double torsion_angle(const std::vector<const PDBAtom*>& atoms,
				    double tol = 0.001);
	
	static double oop_angle(const std::vector<const PDBAtom*>& atoms,
				double tol = 0.001);
	
	static double bending_angle(const std::vector<const PDBAtom*>& atoms);
	
	static double bending_angle(const PDBAtom* atom1,
				    const PDBAtom* atom2,
				    const PDBAtom* atom3);
	
	static double bending_angle(const bmpg_uncc_edu::math::R3Vector<double> & v1,
				    const bmpg_uncc_edu::math::R3Vector<double> & v2,
				    const bmpg_uncc_edu::math::R3Vector<double> & v3);

	static double bond_length(const std::vector<const PDBAtom*>& atoms);
	
	static double torsion_angle_phi(Handle<PDBProteinResidue> res);			//calculates torsion angle phi of the given residue
	static double torsion_angle_psi(Handle<PDBProteinResidue> res);			//calculates torsion angle psi of the given residue
	static double torsion_angle_omega(Handle<PDBProteinResidue> res);		//calculates torsion angle omega of the given residue

	static double torsion_angle(double atom1x,
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
				    double tol = 0.0);

	static double torsion_angle(const bmpg_uncc_edu::math::R3Vector<double> & v1,
				    const bmpg_uncc_edu::math::R3Vector<double> & v2,
				    const bmpg_uncc_edu::math::R3Vector<double> & v3,
				    const bmpg_uncc_edu::math::R3Vector<double> & v4,
				    double tol = 0.0);

	
	static double oop_angle(const bmpg_uncc_edu::math::R3Vector<double> & v1,
				const bmpg_uncc_edu::math::R3Vector<double> & v2,
				const bmpg_uncc_edu::math::R3Vector<double> & v3,
				const bmpg_uncc_edu::math::R3Vector<double> & v4,
				double tol = 0.0);

private:
	static const Element* H_;

};

} // namespace bmpg_uncc_edu::chemistry::helper
} // namespace bmpg_uncc_edu::chemistry
} // namespace bmpg_uncc_edu

#endif

