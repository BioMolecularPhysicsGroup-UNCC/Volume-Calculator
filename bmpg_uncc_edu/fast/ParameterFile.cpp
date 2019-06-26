/**
   Copyright (c) 2008 by Hui Wang
   @author Hui Wang
   @date Jan 20, 2009
*/


#ifndef bmpg_uncc_edu_fast_ParameterFile_cpp
#define bmpg_uncc_edu_fast_ParameterFile_cpp

#include <cstdlib>
#include <algorithm>				//string transform
#include <iostream>
#include <fstream>
#include <sstream>
#include <bmpg_uncc_edu/util/Properties.hpp>
#include <bmpg_uncc_edu/fast/ParameterFile.hpp>


namespace bmpg_uncc_edu {
namespace fast {
using namespace std;
using namespace bmpg_uncc_edu::util;

ParameterFile* ParameterFile::instance_ = NULL;

ParameterFile* ParameterFile::instance()
{
	if(instance_ == NULL){
		instance_ = new ParameterFile();
	}
	return instance_;
}

/**
	The path ends with a "/".
*/
string ParameterFile::path_to_metadata()
{
	string dir;
	if ((dir = props_.property("path_to_metadata")) == Properties::NOT_FOUND) {
		cerr << "\"path_to_metadata\" was not specified in properties file." << endl;
		exit(EXIT_FAILURE);
	}
	if(dir.find_last_of("/") != (dir.size() - 1)){
		dir = dir + "/";
	}
	return dir;
}

string ParameterFile::spectrum_dir()
{
	string dir;
	if ((dir = props_.property("path_to_spectrum")) == Properties::NOT_FOUND) {
		cerr << "\"spectrum_dir\" was not specified in properties file." << endl;
		exit(EXIT_FAILURE);
	}
	if(dir.find_last_of("/") != (dir.size() - 1)){
		dir = dir + "/";
	}
	return dir;
}

string ParameterFile::experimental_data_dir()
{
	string dir;
	if ((dir = props_.property("path_to_experimental_data")) == Properties::NOT_FOUND) {
		cerr << "\"path_to_experimental_data\" was not specified in properties file." << endl;
		exit(EXIT_FAILURE);
	}
	if(dir.find_last_of("/") != (dir.size() - 1)){
		dir = dir + "/";
	}
	return dir;
}

double ParameterFile::energy_bin_size()
{
	string b;
	if ((b = props_.property("energy_bin_size")) == Properties::NOT_FOUND) {
		cerr << "\"energy_bin_size\" was not specified in properties file." << endl;
		exit(EXIT_FAILURE);
	}
	stringstream ss(b);
	double bin_size;
	ss >> bin_size;
	return bin_size;
}


double ParameterFile::T0()
{
	string ts;
	if ((ts = props_.property("ensemble_temperature")) == Properties::NOT_FOUND) {
		cerr << "\"ensemble_temperature\" was not specified in properties file." << endl;
		exit(EXIT_FAILURE);
	}
	stringstream ss(ts);
	double T;
	ss >> T;
	return T;
}

bool ParameterFile::atom_by_occupancy()
{	
	bool atom_by_occupancy = true;
	string s;
	if ((s = props_.property("atom_by_occupancy")) != Properties::NOT_FOUND) {
		std::transform(s.begin(), s.end(), s.begin(), ::tolower);
		atom_by_occupancy = (s == "true");
	}
	return atom_by_occupancy;
}
	


bool ParameterFile::check_for_his()
{
	bool check = true;
	string s;
	if ((s = props_.property("check_for_his")) != Properties::NOT_FOUND) {
		std::transform(s.begin(), s.end(), s.begin(), ::tolower);
		check = (s == "true");
	}
	return check;
}

bool ParameterFile::check_for_css()
{
	bool check = true;
	string s;
	if ((s = props_.property("check_for_css")) != Properties::NOT_FOUND) {
		std::transform(s.begin(), s.end(), s.begin(), ::tolower);
		check = (s == "true");
	}
	return check;
}

bool ParameterFile::calculate_asa()
{
	bool check = true;
	string s;
	if ((s = props_.property("calculate_asa")) != Properties::NOT_FOUND) {
		std::transform(s.begin(), s.end(), s.begin(), ::tolower);
		check = (s == "true");
	}
	return check;
}

bool ParameterFile::calculate_internal_residue_energy()
{
	bool check = true;
	string s;
	if ((s = props_.property("calculate_internal_res_energy")) != Properties::NOT_FOUND) {
		std::transform(s.begin(), s.end(), s.begin(), ::tolower);
		check = (s == "true");
	}
	return check;
}

void ParameterFile::load(const string & config_file_name)
{
	ifstream fin(config_file_name.c_str());
	props_.load(fin, "#!");
	fin.close();
}


}	//namespace bmpg_uncc_edu::fast
}	//namespace bmpg_uncc_edu

#endif

