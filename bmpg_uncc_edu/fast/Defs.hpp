/**
   Copyright (c) 2008 by Hui Wang
   @author Hui Wang, Chuanbin Du
   @date Jun 28, 2009
*/


#ifndef bmpg_uncc_edu_fast_Defs_hpp
#define bmpg_uncc_edu_fast_Defs_hpp

#include <boost/bimap.hpp>
#include <bmpg_uncc_edu/util/Enum2StringBinder.hpp>

namespace bmpg_uncc_edu {
namespace fast {
using namespace bmpg_uncc_edu::util;


/**
	Add alpha and h for hbond and hydrophobic interactions as fudge parameters.
	Add h and sigma for hydration as fudge parameters.
	Add x in accordion idea as a fudge parameter.
	Use Int2Type to collapse types: t1 t2 t3 --> t_s, t4, t5 --> t_c?
*/
typedef enum { TEMPERATURE,
	       PRESSURE,
	       PH,
	       IONIC_CONCENTRATION,
	       COSOLVENT_CONCENTRATION,
	       
	       //solvation states and native state
	       MOBILE,
	       CLATHRATE,
	       BURIED,
	       NATIVE,
	       
	       HBOND,
	       NUM_OF_TORSIONS,
	       NUM_OF_HYDROPHOBIC_CONTACTS,					//A possible secondary order parameter
	       
	       //alpha and h parameters for hbond, hydrophobic, hydration, and vibrational free energy
	       H_HBOND,
	       ALPHA_HBOND,
	       H_HYDROPHOBIC,
	       ALPHA_HYDROPHOBIC,
	       H_HYDRATION,
	       ALPHA_HYDRATION,
	       VIBRATIONAL_REFERENCE_ENERGY,
	       VIBRATIONAL_ENERGY_EXPONENT,
	       
	       NATIVE_SCALE_FACTOR,
	       Q,								//closeness to template structure
	       INDEPENDENT_TORSION_CONSTRAINTS,
//             THETA,                                What is theta (H.Wang) => GLOBAL_FLEXIBILTY  DJJ Dec 8, 2012
	       GLOBAL_FLEXIBILITY,
	       RESIDUE_ENTROPY_SHIFT,
	       INTERACTION_STRAIN_ENERGY,
	       
	       DISORDER_PACKING_ENTROPY,
	       DISORDER_PACKING_DELTA_ENTHALPY,
	       
	       //some constants       
	       UNKNOWN_TYPE
} scanning_variable_t;

typedef enum {
		RIGID_BODY_DOF,
		INTRA_RESIDUE_PASSIVE_BAR,
		HYDRATION_PASSIVE_BAR,
		LINKER_PASSIVE_BAR,
		ZERO_RANK_ENTROPY_CUTOFF,
		ENTROPY_BUNDLING_BINSIZE,	
		LINKER_HYDRATION_RIGIDITY,
		RESIDUE_SOLVENT_ENTROPY_SHIFT,
		RESIDUE_SOLVENT_ENTHALPY_SHIFT,
		
		THREE_STATE_FEL_FLAG,
		NUM_OF_ENERGY_TERMS_IN_MACROSTATE,
		RIGIDITY_SCALING_FACTOR		//the lambda value defined in QSFR notes
} model_constant_t;


bool is_pop(const scanning_variable_t& v);
bool is_entropy_parameter(const scanning_variable_t& v);
bool is_entropy_parameter(const model_constant_t& v);

/**
	The ScanningVariable2String class provides a mapping between strings and the scanning_variable_t.
*/
class ScanningVariable2String : public Enum2StringBinder<scanning_variable_t>
{
public:
	ScanningVariable2String();
	void initialize();

	virtual ~ScanningVariable2String();
};

ostream& operator<<(ostream& os, const scanning_variable_t& rhs);



/**
	The ModelConstant2String class provides a mapping between strings and the model_constant_t.
*/
class ModelConstant2String : public Enum2StringBinder<model_constant_t>
{
public:	
	ModelConstant2String();
	void initialize();

	virtual ~ModelConstant2String();
};

ostream& operator<<(ostream& os, const model_constant_t& rhs);




typedef enum {
	INCLUDE_VIBRATIONAL_ENERGY,
	PLUCKING_FLAG,
	FUDGING_DISORDER_PACKING_INTERACTION,
	SIMPLE_SIGMA_RENORMALIZER,
	INCLUDE_RES_HYDRATION
} model_flag_t;

class ModelFlag2String : public Enum2StringBinder<model_flag_t>
{
public:	
	ModelFlag2String();
	void initialize();

	virtual ~ModelFlag2String();
};

ostream& operator<<(ostream& os, const model_flag_t& rhs);





}	//namespace bmpg_uncc_edu::fast
}	//namespace bmpg_uncc_edu

#endif

