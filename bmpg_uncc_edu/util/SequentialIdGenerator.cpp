/**
   Copyright (c) 2008 by Hui Wang
   @author Hui Wang
   @date Jan 20, 2009
*/


#ifndef bmpg_uncc_edu_fed_graph_traits_details_IDGenerator_cpp
#define bmpg_uncc_edu_fed_graph_traits_details_IDGenerator_cpp

#include <iostream>
#include <bmpg_uncc_edu/util/SequentialIdGenerator.hpp>

namespace bmpg_uncc_edu {
namespace util {
using namespace std;

int SequentialIdGenerator::id_ = 0;

}	//namespace bmpg_uncc_edu::util
}	//namespace bmpg_uncc_edu

#endif

