/**
   Copyright (c) 2008 by Hui Wang
   @author Hui Wang
   @date Mar 12, 2009
*/

#ifndef bmpg_uncc_edu_chemistry_detail_PDBChainManager_cpp
#define bmpg_uncc_edu_chemistry_detail_PDBChainManager_cpp

#include <bmpg_uncc_edu/chemistry/detail/PDBChainManager.hpp>
#include <bmpg_uncc_edu/util/Exception.hpp>

namespace bmpg_uncc_edu {
namespace chemistry {
namespace detail {
using namespace std;

PDBChainManager::PDBChainManager() : curr_chain_id_('\0')
{
}

bool PDBChainManager::has(char c) const
{
	chains_t::const_iterator it;
	it = chains_.find(c);
	return (it == chains_.end());
}

bool PDBChainManager::has(const string& s) const
{
	if(s.empty())
		return false;
	char c = s.at(0);
	return has(c);
}

char PDBChainManager::curr_chain() const{
	return curr_chain_id_;
}

char PDBChainManager::new_chain()
{
	if(curr_chain_id_ == '\0'){
		curr_chain_id_ = 'A';
	} else if(curr_chain_id_ == ' '){
		curr_chain_id_ = 'A';
	} else {
		++curr_chain_id_;
	}
	
	std::pair<chains_t::iterator, bool> result = chains_.insert(curr_chain_id_);
	if(!result.second){
		string s = "Duplicate chain id found; Chain id = " + curr_chain_id_;
		throw bmpg_uncc_edu::util::Exception(s,__FILE__,__LINE__);
	}
	return curr_chain_id_;
}

}	//namespace bmpg_uncc_edu::chemistry::detail
}	//namespace bmpg_uncc_edu::chemistry
}	//namespace bmpg_uncc_edu

#endif

