/**
   Copyright (c) 2008 by Hui Wang
   @author Hui Wang
   @date Jun 03, 2009
*/


#ifndef bmpg_uncc_edu_chemistry_library_loader_AALibLoader_hpp
#define bmpg_uncc_edu_chemistry_library_loader_AALibLoader_hpp

namespace bmpg_uncc_edu {
	namespace chemistry {
		namespace library {
			class AminoAcidLibrary;
		}
	}
}

namespace bmpg_uncc_edu {
namespace chemistry {
namespace library {
namespace loader {

/** 
 * AbstractLoader is an abstract base class for loading libraries.
 */
    class AbstractLoader
{
public:
	virtual void load() = 0;
	virtual ~AbstractLoader(){}
};
/**
 * AALibLoader loads information about amino acids from a file for the AminoAcidLibrary class.
 */
class AALibLoader : public AbstractLoader
{
public:
	AALibLoader(const string& fname,
		    AminoAcidLibrary& lib);
	void load();
private:
	string fname_;
	AminoAcidLibrary& lib_;
};

}	//namespace bmpg_uncc_edu::chemistry::library::loader	
}	//namespace bmpg_uncc_edu::chemistry::library	
}	//namespace bmpg_uncc_edu::chemistry
}	//namespace bmpg_uncc_edu

#endif

