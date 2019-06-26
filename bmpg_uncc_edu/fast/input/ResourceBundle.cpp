/**
   Copyright (c) 2008 by Hui Wang
   @author Hui Wang
   @date Sep 15, 2009
*/


#ifndef bmpg_uncc_edu_fast_input_ResourceBundle_cpp
#define bmpg_uncc_edu_fast_input_ResourceBundle_cpp

#include <iostream>
#include <fstream>
#include <bmpg_uncc_edu/util/Exception.hpp>
#include <bmpg_uncc_edu/fast/input/ResourceBundle.hpp>

namespace bmpg_uncc_edu {
namespace fast {
namespace input {
using namespace std;
using namespace bmpg_uncc_edu::util;

ResourceBundle* ResourceBundle::instance_ = NULL;

ResourceBundle* ResourceBundle::instance()
{
	if(instance_ == NULL){
		instance_ = new ResourceBundle();
	}
	return instance_;
}

ResourceBundle::ResourceBundle()
{
	
}
/*
void ResourceBundle::load(const string& fname)
{
	ifstream file;
	file.open(fname.c_str());
	Properties::load(file);
}
*/
string ResourceBundle::get(const string& key) const
{
	string result = property(key);
	if(result.compare(bmpg_uncc_edu::util::Properties::NOT_FOUND) == 0){
		throw Exception("The key " + key + " is not found in the resources.", __FILE__, __LINE__);	
	}
	return result;
}

}	//namespace bmpg_uncc_edu::fast::input
}	//namespace bmpg_uncc_edu::fast
}	//namespace bmpg_uncc_edu

#endif

