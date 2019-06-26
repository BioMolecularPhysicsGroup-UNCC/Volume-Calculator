/**
   Copyright (c) 2008 by Hui Wang
   @author Hui Wang, Chuanbin Du
   @date Jun 28, 2009
*/


#ifndef bmpg_uncc_edu_fast_Defs_cpp
#define bmpg_uncc_edu_fast_Defs_cpp

#include <sstream>
#include <bmpg_uncc_edu/fast/Defs.hpp>
#include <bmpg_uncc_edu/fast/input/ResourceBundle.hpp>
#include <bmpg_uncc_edu/util/Exception.hpp>


namespace bmpg_uncc_edu {
namespace fast {
using namespace std;
using namespace bmpg_uncc_edu::fast::input;
using namespace bmpg_uncc_edu::util;

ScanningVariable2String::ScanningVariable2String()
{
}

ScanningVariable2String::~ScanningVariable2String()
{
}

void ScanningVariable2String::initialize()
{
	ResourceBundle* resources = ResourceBundle::instance();
	
	insert(value_t(TEMPERATURE,resources->property("TEMPERATURE")));
	insert(value_t(PRESSURE,resources->property("PRESSURE")));
	insert(value_t(PH,resources->property("PH")));
	insert(value_t(IONIC_CONCENTRATION,resources->property("IONIC_CONCENTRATION")));
	insert(value_t(COSOLVENT_CONCENTRATION,resources->property("COSOLVENT_CONCENTRATION")));
	insert(value_t(MOBILE,resources->property("MOBILE")));	
	insert(value_t(CLATHRATE,resources->property("CLATHRATE")));
	insert(value_t(BURIED,resources->property("BURIED")));
	insert(value_t(NATIVE,resources->property("NATIVE")));
	insert(value_t(HBOND,resources->property("HBOND")));
	insert(value_t(NUM_OF_TORSIONS,resources->property("NUM_OF_TORSIONS")));
	insert(value_t(NUM_OF_HYDROPHOBIC_CONTACTS,resources->property("NUM_OF_HYDROPHOBIC_CONTACTS")));
	
	
	
	insert(value_t(H_HBOND,resources->property("H_HBOND")));
	insert(value_t(ALPHA_HBOND,resources->property("ALPHA_HBOND")));
	insert(value_t(H_HYDROPHOBIC,resources->property("H_HYDROPHOBIC")));
	insert(value_t(ALPHA_HYDROPHOBIC,resources->property("ALPHA_HYDROPHOBIC")));
	insert(value_t(H_HYDRATION,resources->property("H_HYDRATION")));

	insert(value_t(VIBRATIONAL_REFERENCE_ENERGY,resources->property("VIBRATIONAL_REFERENCE_ENERGY")));
	insert(value_t(VIBRATIONAL_ENERGY_EXPONENT,resources->property("VIBRATIONAL_ENERGY_EXPONENT")));


	insert(value_t(ALPHA_HYDRATION,resources->property("ALPHA_HYDRATION")));
	insert(value_t(NATIVE_SCALE_FACTOR,resources->property("NATIVE_SCALE_FACTOR")));
	
	insert(value_t(Q,resources->property("Q")));
	insert(value_t(INDEPENDENT_TORSION_CONSTRAINTS,resources->property("INDEPENDENT_TORSION_CONSTRAINTS")));
//	insert(value_t(THETA,resources->property("THETA")));     to be removed completely: DJJ Dec 8, 2012
	insert(value_t(GLOBAL_FLEXIBILITY,resources->property("GLOBAL_FLEXIBILITY")));
	
	
	insert(value_t(RESIDUE_ENTROPY_SHIFT,resources->property("RESIDUE_ENTROPY_SHIFT")));
	
	insert(value_t(INTERACTION_STRAIN_ENERGY,resources->property("INTERACTION_STRAIN_ENERGY")));
	
	insert(value_t(DISORDER_PACKING_ENTROPY,resources->property("DISORDER_PACKING_ENTROPY")));
	insert(value_t(DISORDER_PACKING_DELTA_ENTHALPY,resources->property("DISORDER_PACKING_DELTA_ENTHALPY")));
	insert(value_t(DISORDER_PACKING_DELTA_ENTHALPY,resources->property("DISORDER_PACKING_DELTA_ENTHALPY")));

}

bool is_pop(const scanning_variable_t& v)
{
	return (v == BURIED) || (v == MOBILE) || (v == CLATHRATE) || (v == NATIVE);
}

bool is_entropy_parameter(const scanning_variable_t& v)
{
	return (v == ALPHA_HBOND || v == ALPHA_HYDROPHOBIC || 
		v == ALPHA_HYDRATION || v == RESIDUE_ENTROPY_SHIFT ||
		v == DISORDER_PACKING_ENTROPY);
}

bool is_entropy_parameter(const model_constant_t& v)
{
	return (v == ZERO_RANK_ENTROPY_CUTOFF || v == ENTROPY_BUNDLING_BINSIZE || v == RESIDUE_SOLVENT_ENTROPY_SHIFT);
}

		
ostream& operator<<(ostream& os, const scanning_variable_t& rhs)
{
	ResourceBundle* resources = ResourceBundle::instance();
	
	string s;
	switch(rhs){
		case TEMPERATURE:
			s = resources->property("TEMPERATURE");
			break;
		case PRESSURE:
			s = resources->property("PRESSURE");
			break;
		case PH:
			s = resources->property("PH");
			break;
		case IONIC_CONCENTRATION:
			s = resources->property("IONIC_CONCENTRATION");
			break;
		case COSOLVENT_CONCENTRATION:
			s = resources->property("COSOLVENT_CONCENTRATION");
			break;
		case MOBILE:
			s = resources->property("MOBILE");
			break;
		case CLATHRATE:
			s = resources->property("CLATHRATE");
			break;
		case BURIED:
			s = resources->property("BURIED");
			break;
		case NATIVE:
			s = resources->property("NATIVE");
			break;
		case HBOND:
			s = resources->property("HBOND");
			break;
		case NUM_OF_TORSIONS:
			s = resources->property("NUM_OF_TORSIONS");
			break;
		case NUM_OF_HYDROPHOBIC_CONTACTS:
			s = resources->property("NUM_OF_HYDROPHOBIC_CONTACTS");
			break;
		case Q:
			s = resources->property("Q");
			break;
		case INDEPENDENT_TORSION_CONSTRAINTS:
			s = resources->property("INDEPENDENT_TORSION_CONSTRAINTS");
			break;
/*              to be removed completely: DJJ Dec 8, 2012 
		case THETA:
			s = resources->property("THETA");
			break;
*/
		case H_HBOND:
			s = resources->property("H_HBOND");
			break;
		case ALPHA_HBOND:
			s = resources->property("ALPHA_HBOND");
			break;
		case H_HYDROPHOBIC:
			s = resources->property("H_HYDROPHOBIC");
			break;
		case ALPHA_HYDROPHOBIC:
			s = resources->property("ALPHA_HYDROPHOBIC");
			break;
		case H_HYDRATION:
			s = resources->property("H_HYDRATION");
			break;
		case ALPHA_HYDRATION:
			s = resources->property("ALPHA_HYDRATION");
			break;
		case VIBRATIONAL_REFERENCE_ENERGY:
			s = resources->property("VIBRATIONAL_REFERENCE_ENERGY");
			break;
		case VIBRATIONAL_ENERGY_EXPONENT:
			s = resources->property("VIBRATIONAL_ENERGY_EXPONENT");
			break;
		case NATIVE_SCALE_FACTOR:
			s = resources->property("NATIVE_SCALE_FACTOR");
			break;
		case GLOBAL_FLEXIBILITY:
			s = resources->property("GLOBAL_FLEXIBILITY");
			break;
		case RESIDUE_ENTROPY_SHIFT:
			s = resources->property("RESIDUE_ENTROPY_SHIFT");
			break;	
		case INTERACTION_STRAIN_ENERGY:
			s = resources->property("INTERACTION_STRAIN_ENERGY");
			break;	
		case DISORDER_PACKING_ENTROPY:
			s = resources->property("DISORDER_PACKING_ENTROPY");
			break;	
		case DISORDER_PACKING_DELTA_ENTHALPY:
			s = resources->property("DISORDER_PACKING_DELTA_ENTHALPY");
			break;		
		default:
			s = "UNRECOGNIZED_TYPE";
			break;
	}

	os << s;
	return os;	
}



ModelConstant2String::ModelConstant2String()
{
}

ModelConstant2String::~ModelConstant2String()
{
}

void ModelConstant2String::initialize()
{
	ResourceBundle* resources = ResourceBundle::instance();
	
	insert(value_t(RIGID_BODY_DOF,resources->property("RIGID_BODY_DOF")));
	insert(value_t(INTRA_RESIDUE_PASSIVE_BAR,resources->property("INTRA_RESIDUE_PASSIVE_BAR")));
	insert(value_t(HYDRATION_PASSIVE_BAR,resources->property("HYDRATION_PASSIVE_BAR")));
	insert(value_t(LINKER_PASSIVE_BAR,resources->property("LINKER_PASSIVE_BAR")));
	insert(value_t(ZERO_RANK_ENTROPY_CUTOFF,resources->property("ZERO_RANK_ENTROPY_CUTOFF")));
	insert(value_t(ENTROPY_BUNDLING_BINSIZE,resources->property("ENTROPY_BUNDLING_BINSIZE")));
	insert(value_t(LINKER_HYDRATION_RIGIDITY,resources->property("LINKER_HYDRATION_RIGIDITY")));
	
	insert(value_t(RESIDUE_SOLVENT_ENTROPY_SHIFT,resources->property("RESIDUE_SOLVENT_ENTROPY_SHIFT")));
	insert(value_t(RESIDUE_SOLVENT_ENTHALPY_SHIFT,resources->property("RESIDUE_SOLVENT_ENTHALPY_SHIFT")));
		
	insert(value_t(THREE_STATE_FEL_FLAG,resources->property("THREE_STATE_FEL_FLAG")));	
	
	insert(value_t(NUM_OF_ENERGY_TERMS_IN_MACROSTATE,resources->property("NUM_OF_ENERGY_TERMS_IN_MACROSTATE")));
	insert(value_t(RIGIDITY_SCALING_FACTOR,resources->property("RIGIDITY_SCALING_FACTOR")));
}


ostream& operator<<(ostream& os, const model_constant_t& rhs)
{
	ResourceBundle* resources = ResourceBundle::instance();
	
	string s;
	switch(rhs){
		case RIGID_BODY_DOF:
			s = resources->property("RIGID_BODY_DOF");
			break;
		case INTRA_RESIDUE_PASSIVE_BAR:
			s = resources->property("INTRA_RESIDUE_PASSIVE_BAR");
			break;
		case HYDRATION_PASSIVE_BAR:
			s = resources->property("HYDRATION_PASSIVE_BAR");
			break;
		case LINKER_PASSIVE_BAR:
			s = resources->property("LINKER_PASSIVE_BAR");
			break;
		case ZERO_RANK_ENTROPY_CUTOFF:
			s = resources->property("ZERO_RANK_ENTROPY_CUTOFF");
			break;
		case ENTROPY_BUNDLING_BINSIZE:
			s = resources->property("ENTROPY_BUNDLING_BINSIZE");
			break;	
		case THREE_STATE_FEL_FLAG:
			s = resources->property("THREE_STATE_FEL_FLAG");
			break;	
		case NUM_OF_ENERGY_TERMS_IN_MACROSTATE:
			s = resources->property("NUM_OF_ENERGY_TERMS_IN_MACROSTATE");
			break;
		case RIGIDITY_SCALING_FACTOR:
			s = resources->property("RIGIDITY_SCALING_FACTOR");
			break;
		case LINKER_HYDRATION_RIGIDITY:
			s = resources->property("LINKER_HYDRATION_RIGIDITY");
			break;	
		case RESIDUE_SOLVENT_ENTROPY_SHIFT:
			s = resources->property("RESIDUE_SOLVENT_ENTROPY_SHIFT");
			break;
		case RESIDUE_SOLVENT_ENTHALPY_SHIFT:
			s = resources->property("RESIDUE_SOLVENT_ENTHALPY_SHIFT");
			break;	
		default:
			s = "UNRECOGNIZED_TYPE";
			break;
	}

	os << s;
	return os;	
}



ModelFlag2String::ModelFlag2String()
{
}

ModelFlag2String::~ModelFlag2String()
{
}

void ModelFlag2String::initialize()
{
	ResourceBundle* resources = ResourceBundle::instance();
	insert(value_t(INCLUDE_VIBRATIONAL_ENERGY,resources->property("INCLUDE_VIBRATIONAL_ENERGY")));
	insert(value_t(PLUCKING_FLAG,resources->property("PLUCKING_FLAG")));	
	insert(value_t(FUDGING_DISORDER_PACKING_INTERACTION,resources->property("FUDGING_DISORDER_PACKING_INTERACTION")));
	insert(value_t(SIMPLE_SIGMA_RENORMALIZER,resources->property("SIMPLE_SIGMA_RENORMALIZER")));
	insert(value_t(INCLUDE_RES_HYDRATION,resources->property("INCLUDE_RES_HYDRATION")));
}


ostream& operator<<(ostream& os, const model_flag_t& rhs)
{
	ResourceBundle* resources = ResourceBundle::instance();
	
	string s;
	switch(rhs){
		case INCLUDE_VIBRATIONAL_ENERGY:
			s = resources->property("INCLUDE_VIBRATIONAL_ENERGY");
			break;
		case PLUCKING_FLAG:
			s = resources->property("PLUCKING_FLAG");
			break;
			
		case SIMPLE_SIGMA_RENORMALIZER:
			s = resources->property("SIMPLE_SIGMA_RENORMALIZER");
			break;
		case FUDGING_DISORDER_PACKING_INTERACTION:
			s = resources->property("FUDGING_DISORDER_PACKING_INTERACTION");
			break;
			
		case INCLUDE_RES_HYDRATION:
			s = resources->property("INCLUDE_RES_HYDRATION");
			break;
			
		default:
			s = "UNRECOGNIZED_TYPE";
			break;
	}
	return os;
}

}	//namespace bmpg_uncc_edu::fast
}	//namespace bmpg_uncc_edu

#endif

