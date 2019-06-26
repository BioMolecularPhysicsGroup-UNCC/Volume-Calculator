/**
   Copyright (c) 2008 by Hui Wang
   @author Hui Wang
   @date Jan 20, 2009
*/


#ifndef bmpg_uncc_edu_fast_ParameterFile_hpp
#define bmpg_uncc_edu_fast_ParameterFile_hpp

#include <string>
#include <bmpg_uncc_edu/util/Properties.hpp>

namespace bmpg_uncc_edu {
namespace fast {
using namespace std;
using namespace bmpg_uncc_edu::util;

/**
 * ParameterFile loads the properties configuration file and provides wrapper methods to access the information stored there.
 * This class is specific to a particular implementation of the DCM and can be modified as necessary to accommodate different models.
 */
class ParameterFile
{
public:
	static ParameterFile* instance();
	void load(const string & config_file_name);
	
	double T0();
	double energy_bin_size();

	bool calculate_internal_residue_energy();
	bool calculate_asa();
//	bool calculate_inter_residue_energy();           not implemented
	
	bool check_for_his();
	bool check_for_css();

	bool atom_by_occupancy();
	
	string spectrum_dir();
	string path_to_metadata();
	string experimental_data_dir();

private:
	static ParameterFile* instance_;
	Properties props_;
private:
	ParameterFile(){}
	ParameterFile(const ParameterFile&);
	ParameterFile& operator=(const ParameterFile&);	
};

}	//namespace bmpg_uncc_edu::fast
}	//namespace bmpg_uncc_edu

#endif

