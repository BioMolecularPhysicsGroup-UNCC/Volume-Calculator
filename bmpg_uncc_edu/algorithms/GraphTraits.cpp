/**
   Copyright (c) 2008 by Hui Wang
   @author Hui Wang
   @date Jan 20, 2009
*/


#ifndef bmpg_uncc_edu_algorithms_Graph_Traits_cpp
#define bmpg_uncc_edu_algorithms_Graph_Traits_cpp

#include <iostream>
#include <string>
#include <limits>
#include <list>
#include <map>
#include <boost/shared_ptr.hpp>
#include <bmpg_uncc_edu/algorithms/GraphTraits.hpp>


namespace bmpg_uncc_edu {
namespace algorithms {
using namespace std;



int Limit<int>::inf()
{
	return numeric_limits<int>::max();	
}


double Limit<double>::inf()
{
	return numeric_limits<double>::max();
}


float Limit<float>::inf()
{
	return numeric_limits<float>::max();	
}

int Limit<int>::zero()
{
	return 0;
}

double Limit<double>::zero()
{
	return 0;	
}

float Limit<float>::zero()
{
	return 0;
}

}	//namespace bmpg_uncc_edu::algorithms
}	//namespace bmpg_uncc_edu

#endif

