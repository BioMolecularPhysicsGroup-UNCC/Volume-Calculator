/**
   Copyright (c) 2008 by Mike Fairchild, Hui Wang
   @author: Mike Fairchild, http://www.mikef.org, mjfairch@gmail.com
   @author: Hui Wang
*/

#ifndef bmpg_uncc_edu_chemistry_Atom_hpp
#define bmpg_uncc_edu_chemistry_Atom_hpp

#include <string>

namespace bmpg_uncc_edu {
namespace chemistry {
using namespace std;

class Element;

/**
 * Atom represents a single atom in the periodic table, along with the mass, and the x, y, and z coordinates of the atom.
 */

class Atom {
public:
	double x;								// X coordinate in Angstroms
	double y;								// Y coordinate in Angstroms
	double z;								// Z coordinate in Angstroms
	
	Atom();
	Atom(const Atom& rhs);

	explicit Atom(const string & name);
	Atom(const Element * elem,
	     double x,
	     double y,
	     double z);
	
	size_t Z() const;
	
	virtual double mass() const;
	
	const Element * element() const { return element_; }
	
	void element(const Element * e);
	
	Atom& operator=(const Atom &rhs);
	
	friend ostream& operator<<(ostream & out, const Atom &rhs);
	
	virtual ~Atom(){}
protected:
	const Element * element_;
};

} // namespace bmpg_uncc_edu::chemistry
} // namespace bmpg_uncc_edu

#endif

