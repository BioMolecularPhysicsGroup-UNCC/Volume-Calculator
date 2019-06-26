/**
  Copyright (c) 2008 by Mike Fairchild
  @author: Mike Fairchild, http://www.mikef.org, mjfairch@gmail.com
*/

#ifndef bmpg_uncc_edu_chemistry_library_PeriodicTable_hpp
#define bmpg_uncc_edu_chemistry_library_PeriodicTable_hpp

#include <map>

namespace bmpg_uncc_edu {
	namespace chemistry {
		class Element;
	}
}

namespace bmpg_uncc_edu {
namespace chemistry {
namespace library {
	

/**
 * PeriodicTable stores the collection of Element classes found in the periodic table.
 */    

class PeriodicTable {
	public:
		typedef enum { H=1, He, Li, Be, B, C, N, O, F, Ne, Na, Mg, Al, Si,
			P, S, Cl, Ar, K, Ca, Sc, Ti, V, Cr, Mn, Fe, Co, Ni, Cu, Zn,
			Ga, Ge, As, Se, Br, Kr, Rb, Sr, Y, Zr, Nb, Mo, Tc, Ru, Rh,
			Pd, Ag, Cd, In, Sn, Sb, Te, I, Xe, Cs, Ba, La, Ce, Pr, Nd,
			Pm, Sm, Eu, Gd, Tb, Dy, Ho, Er, Tm, Yb, Lu, Hf, Ta, W, Re,
			Os, Ir, Pt, Au, Hg, Tl, Pb, Bi, Po, At, Rn, Fr, Ra, Ac, Th,
			Pa, U, Np, Pu, Am, Cm, Bk, Cf, Es, Fm, Md, No, Lr, Rf, Db,
			Sg, Bh, Hs, Mt, _110 } atomic_symbol_t;

		// Periodic Table is a singleton.  This method provides access to the
		// one and only copy of it.
		static PeriodicTable * instance();
		// Returns a pointer to the most abundant element with atomic number Z.
		// If A isn't specified or is zero, then the most abundant isotope with
		// atomic number Z is returned.
		const Element * const element(size_t Z, size_t A = 0) const;
		
		const Element * const element(const std::string & symbol,
					      size_t A = 0) const;
		
		std::istream & load(std::istream & is); // Loads the periodic table data from the given input stream		
		friend std::istream & operator>>(std::istream & is, PeriodicTable & pt);
		friend std::ostream & operator<<(std::ostream & os, const PeriodicTable & pt);

	private:
		// Hide constructor, copy constructor, operator=, and destructor
		PeriodicTable();
		PeriodicTable(const PeriodicTable & rhs);
		PeriodicTable & operator=(const PeriodicTable & rhs);
		~PeriodicTable();
		
	private:
		static PeriodicTable * instance_;
		typedef std::map<std::pair<size_t, size_t>, Element *> table_t;
		table_t table_; // Maps (Z,A) to Element *
};

}	//namespace bmpg_uncc_edu::chemistry::library
} 	// namespace bmpg_uncc_edu::chemistry
} 	// namespace bmpg_uncc_edu

#endif

