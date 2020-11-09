/**
   Copyright (c) 2008 by Hui Wang
   Modified by Jenny Farmer 2019
   @author Hui Wang
   @date Jan 20, 2009
*/


#ifndef bmpg_uncc_edu_fast_FASTInitializer_hpp
#define bmpg_uncc_edu_fast_FASTInitializer_hpp

#include <string>

using namespace std;

/** 
 * FASTInitializer loads all the user preferences for FAST.  The configuration file specifies the user preferences, as well as paths and
 * file names necessary to load all required libraries. Basic error checking is performed to ensure FAST has minimal information to run.
 */

class FASTInitializer
{
public:
	FASTInitializer();
	void load_user_defined_libraries();
	void close();
};


#endif

