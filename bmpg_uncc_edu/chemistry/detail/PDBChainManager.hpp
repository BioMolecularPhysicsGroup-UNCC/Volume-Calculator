/**
   Copyright (c) 2008 by Hui Wang
   @author Hui Wang
   @date Mar 12, 2009
*/

#ifndef bmpg_uncc_edu_chemistry_detail_PDBChainManager_hpp
#define bmpg_uncc_edu_chemistry_detail_PDBChainManager_hpp

#include <set>
#include <string>
#include <cstring>

namespace bmpg_uncc_edu {
namespace chemistry {
namespace detail {
using namespace std;

class PDBChainManager
{
public:
	typedef std::set<char> chains_t;
	PDBChainManager();
	bool has(char c) const;
	bool has(const string& s) const;
	char new_chain();
	char curr_chain() const;

private:
	char curr_chain_id_;
	chains_t chains_;
};

}	//namespace bmpg_uncc_edu::chemistry::detail
}	//namespace bmpg_uncc_edu::chemistry
}	//namespace bmpg_uncc_edu

#endif

