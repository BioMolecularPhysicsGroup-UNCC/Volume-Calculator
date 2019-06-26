/**
   Copyright (c) 2008 by Hui Wang
   @author Hui Wang
   @date Jan 20, 2009
*/


#ifndef bmpg_uncc_edu_chemistry_CovalentBond_cpp
#define bmpg_uncc_edu_chemistry_CovalentBond_cpp

#include <iostream>
#include <bmpg_uncc_edu/chemistry/CovalentBond.hpp>
#include <bmpg_uncc_edu/chemistry/Bond.hpp>
#include <bmpg_uncc_edu/chemistry/PDBAtom.hpp>

namespace bmpg_uncc_edu {
namespace chemistry {
using namespace std;

CovalentBond::CovalentBond(const PDBAtom* atom1,
			   const PDBAtom* atom2) :
			   first_atom_(atom1),
			   second_atom_(atom2)
{
}

const PDBAtom* CovalentBond::first_atom() const
{
	return first_atom_;
}
	
const PDBAtom* CovalentBond::second_atom() const
{
	return second_atom_;
}

ostream& operator<<(ostream & out, const CovalentBond & rhs)
{
	out << *rhs.first_atom() << "--" << *rhs.second_atom();
	return out;
}

}	//namespace bmpg_uncc_edu::chemistry
}	//namespace bmpg_uncc_edu

#endif

