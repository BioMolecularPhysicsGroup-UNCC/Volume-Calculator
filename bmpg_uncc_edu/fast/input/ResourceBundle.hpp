/**
   Copyright (c) 2008 by Hui Wang
   @author Hui Wang
   @date Sep 15, 2009
*/


#ifndef bmpg_uncc_edu_fast_input_ResourceBundle_hpp
#define bmpg_uncc_edu_fast_input_ResourceBundle_hpp

#include <string>
#include <bmpg_uncc_edu/util/Properties.hpp>

namespace bmpg_uncc_edu {
namespace fast {
namespace input {
using namespace std;
using namespace bmpg_uncc_edu::util;

/**
 * ResourceBundle inherits from the bmpg_uncc_edu::util::Properties file. It is intended to
 * provide a mapping between the keywords in our program and the strings that are used in the input file.
*/
class ResourceBundle : public Properties
{
public:	
	ResourceBundle();
	
	string get(const string& key) const;
	static ResourceBundle* instance();
private:
	static ResourceBundle* instance_;
	ResourceBundle(const ResourceBundle& rhs);
	ResourceBundle& operator=(const ResourceBundle& rhs);
};

}	//namespace bmpg_uncc_edu::fast::input
}	//namespace bmpg_uncc_edu::fast
}	//namespace bmpg_uncc_edu

#endif

