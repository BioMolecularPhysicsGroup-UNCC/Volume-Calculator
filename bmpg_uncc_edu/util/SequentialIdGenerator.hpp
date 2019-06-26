/**
   Copyright (c) 2008 by Hui Wang
   @author Hui Wang
   @date Mar 19, 2010
*/


#ifndef bmpg_uncc_edu_util_SequentialIdGenerator_hpp
#define bmpg_uncc_edu_util_SequentialIdGenerator_hpp

#include <cstdlib>

namespace bmpg_uncc_edu {
namespace util {	
using namespace std;

/**
	Generate integer id's sequentially, starting from 0.
*/
class SequentialIdGenerator
{
public:
	SequentialIdGenerator() {}
	static int get() {return id_++;}
private:
	static int id_;
};


}	//namespace bmpg_uncc_edu::util
}	//namespace bmpg_uncc_edu

#endif

