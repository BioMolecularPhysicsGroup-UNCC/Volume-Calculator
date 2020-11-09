/* 
 * File:   AtomCoordinates.hpp
 * Author: jenny
 *
 * Created on November 18, 2016, 12:50 PM
 */

#ifndef ATOMCOORDINATES_HPP
#define	ATOMCOORDINATES_HPP

#include <boost/random.hpp>
#include <fstream>

#include <bmpg_uncc_edu/fast/FASTInitializer.hpp>
//#include "FASTInitializer.hpp";

#include <bmpg_uncc_edu/util/Exception.hpp>
//#include "PDBProtein.hpp"
#include "bmpg_uncc_edu/chemistry/PDBProtein.hpp"
#include "bmpg_uncc_edu/chemistry/PDBAtom.hpp"
#include "bmpg_uncc_edu/chemistry/Element.hpp"
#include <bmpg_uncc_edu/util/logger/LoggerFactory.hpp>
#include <bmpg_uncc_edu/chemistry/helper/PDBHelper.hpp>

using namespace std;
using namespace bmpg_uncc_edu::util;
using namespace bmpg_uncc_edu::chemistry;
using namespace bmpg_uncc_edu::chemistry::hbond;
using namespace bmpg_uncc_edu::chemistry::helper;
using namespace bmpg_uncc_edu::util::logger;
using namespace bmpg_uncc_edu::algorithms;

class AtomCoordinates {
public:
    AtomCoordinates(bool bondi);
    virtual ~AtomCoordinates();
    
    void getAtomCoordinates(string PDBFileName);
    void rotateCoordinates();
    void residueVolumeSum(int * cavityPerAtom, int * microvoidPerAtom, int * boundaryPerAtom, int * proteinPerAtom, int * cavityPerAtomClose, int * mvPerAtomClose, int * boundaryPerAtomClose);
    void writeFractionalVolumes(int nRotations, const char * fileNameResidue, const char * fileNameAtom, const char * fileNameResidueClose, const char * fileNameAtomClose);
      
    int    nAtoms;      
    int    nResidues;
    float ** atomCoords;
    
private:
    bool    bondi_;        
    float hydrogenVDW_ = 0;
    
    const char  ** atomTypes;
    const char  ** residueTypes;
    
    char * chains_;
    int    nChains_;
    int  * maxChainResidues_;
    int    maxResiduesInChain_;
    const int   maxChainsAllowed_ = 8;
        
    float  ** cavityPerResidue1_;
    float  ** microvoidPerResidue1_;
    float  ** boundaryPerResidue1_;
    float  ** proteinPerResidue1_;
    
         
    float  ** cavityPerResidueClose1_;
    float  ** mvPerResidueClose1_;
    float  ** boundaryPerResidueClose1_;

    float  ** cavityPerResidue2_;                                                //used as second moments
    float  ** microvoidPerResidue2_;
    float  ** boundaryPerResidue2_;
    float  ** proteinPerResidue2_;
    
    float  ** cavityPerResidueClose2_;    
    float  ** mvPerResidueClose2_;
    float  ** boundaryPerResidueClose2_;
    
    float  * cavityPerAtom1_;
    float  * microvoidPerAtom1_;
    float  * boundaryPerAtom1_;
    float  * proteinPerAtom1_;
    
    float  * cavityPerAtomClose1_;
    float  * mvPerAtomClose1_;
    float  * boundaryPerAtomClose1_;

    float  * cavityPerAtom2_;                                               
    float  * microvoidPerAtom2_;
    float  * boundaryPerAtom2_;
    float  * proteinPerAtom2_;
    
    float  * cavityPerAtomClose2_;     
    float  * mvPerAtomClose2_; 
    float  * boundaryPerAtomClose2_;
    
    float rotatex(float phi, float theta, float psi, float x, float y, float z);
    float rotatey(float phi, float theta, float psi, float x, float y, float z);
    float rotatez(float phi, float theta, float psi, float x, float y, float z);
    void alignAxes();
    int dsyevj3(float A[3][3], float Q[3][3], float w[3]);
};

#endif	

