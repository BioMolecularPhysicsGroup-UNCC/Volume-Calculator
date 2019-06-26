/**
   Copyright (c) 2008 by Mike Fairchild, Hui Wang
   @author: Mike Fairchild, http://www.mikef.org, mjfairch@gmail.com
   @author: Hui Wang
*/

#ifndef bmpg_uncc_edu_chemistry_Atom_cpp
#define bmpg_uncc_edu_chemistry_Atom_cpp

#include <iostream>
#include <bmpg_uncc_edu/chemistry/Element.hpp>
#include <bmpg_uncc_edu/chemistry/Atom.hpp>
#include <bmpg_uncc_edu/chemistry/library/PeriodicTable.hpp>

namespace bmpg_uncc_edu {
namespace chemistry {
using namespace std;

Atom::Atom() : element_(NULL)
{
	x = y = z = 0;
}
Atom::Atom(const string & name)
{
	library::PeriodicTable *pt = library::PeriodicTable::instance();
	element_ = pt->element(name);
	x = y = z = 0;
}
Atom::Atom(const Element * elem, double xpos, double ypos, double zpos) : x(xpos), y(ypos), z(zpos), element_(elem)
{

}

Atom::Atom(const Atom& rhs)
{
	operator=(rhs);
}

size_t Atom::Z() const
{
	return element_->Z();
}
	
double Atom::mass() const
{
	return element_->atomic_mass();
}	

void Atom::element(const Element * e)
{
	element_ = e;
}

Atom& Atom::operator=(const Atom &rhs)
{
	element_ = rhs.element_;
	return *this;
}

	
ostream& operator<<(ostream & out, const Atom &rhs)
{
	out << rhs.element_->Z();
	return out;
}

} // namespace bmpg_uncc_edu::chemistry
} // namespace bmpg_uncc_edu

#endif

