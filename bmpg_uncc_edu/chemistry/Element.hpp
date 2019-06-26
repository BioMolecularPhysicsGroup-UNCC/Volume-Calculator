/**
   Copyright (c) 2008 by Mike Fairchild, Hui Wang
   @author: Mike Fairchild, http://www.mikef.org, mjfairch@gmail.com
   @author: Hui Wang
*/

#ifndef bmpg_uncc_edu_chemistry_Element_hpp
#define bmpg_uncc_edu_chemistry_Element_hpp

#include <cstddef>
#include <string>
#include <vector>

namespace bmpg_uncc_edu {
namespace chemistry {


/**
 * Element contains all the information about a single element in the periodic table, 
 * such as name, atomic number, atomic mass, half-life, etc.
 */    

class Element
{
	friend class PeriodicTable;
	
	public:
		Element() : Z_(0), A_(0), N_(0), symbol_(""), name_(""), 
			chemical_weight_(0.0), atomic_mass_(0.0), abundance_(0.0), 
			radioactive_(false), half_life_(0.0)
			{
			}
		// default copy constructor, operator=, and destructor

		enum DECAY_MODE { UNKNOWN_DECAY=1, ALPHA_DECAY, BETA_MINUS_DECAY,
			BETA_PLUS_DECAY, ELECTRON_CAPTURE_DECAY };
		
		size_t Z() const;
		inline size_t A() const { return A_; }
		inline size_t N() const { return N_; }
		inline std::string symbol() const { return symbol_; }
		inline std::string name() const { return name_; }
		inline double chemical_weight() const { return chemical_weight_; }
		inline double atomic_mass() const { return atomic_mass_; }
		inline double abundance() const { return abundance_; }
                                                                                // Added electronegativity 1/21/2017
                inline double electronegativity() const { return electronegativity_; }  
  		inline bool radioactive() const { return radioactive_; }
		inline double half_life() const { return half_life_; }
		inline const std::vector<DECAY_MODE> & decay_modes() const { return decay_modes_; }
		
		friend std::ostream & operator<<(std::ostream & os, const Element & e);
		friend std::istream & operator>>(std::istream & is, Element & e);

	private:
		size_t Z_; 							// Atomic number (number of protons)
		size_t A_; 							// Mass number (number of protons and neutrons)
		size_t N_; 							// Number of neutrons (equals A - Z)
		std::string symbol_; 						// atomic symbol
		std::string name_; 						// element name
		double chemical_weight_; 					// Weighted (average) mass (in u) of element across different isotopes
		double atomic_mass_; 						// Atomic mass of this isotope
		double abundance_; 						// Fractional abundance of this isotope found naturally
                double electronegativity_;                                      // Added electronegativity 1/21/2017
		bool radioactive_; 						// Is this isotope radioactive?
		double half_life_;						// Half-life of this isotope in seconds (s)
		std::vector<DECAY_MODE> decay_modes_; 				// A list of the decay modes through which this isotope decays
};

} // namespace bmpg_uncc_edu::chemistry
} // namespace bmpg_uncc_edu

#endif

